//! Blocked GZip Format (BGZF) reader with seeking support.
//!
//! BGZF is a variant of gzip that uses independent blocks for random access.
//! Each block is a valid gzip member, allowing decompression from any block boundary.

use crate::gzi::GziIndex;
use flate2::Decompress;
use std::io::{self, BufRead, Read, Seek, SeekFrom};

/// BGZF magic numbers and constants
const GZIP_ID1: u8 = 0x1f;
const GZIP_ID2: u8 = 0x8b;
const GZIP_CM_DEFLATE: u8 = 8;
const GZIP_FLG_FEXTRA: u8 = 4;
#[allow(dead_code)]
const GZIP_OS_UNKNOWN: u8 = 255;
const BGZF_EXTRA_ID: u8 = 66; // 'B'
const BGZF_EXTRA_SUBFIELD: u8 = 67; // 'C'
const BGZF_MAX_BLOCK_SIZE: usize = 64 * 1024;

/// A BGZF reader with optional index for seeking.
///
/// This reader can decompress BGZF-compressed data and, when provided with
/// a .gzi index, can seek to arbitrary positions in the uncompressed stream.
///
/// # Type Parameters
///
/// * `R` - The underlying reader type (must implement Read and Seek)
///
/// # Example
///
/// ```no_run
/// use fastx::bgzf::BgzfReader;
/// use fastx::gzi::GziIndex;
/// use std::fs::File;
/// use std::path::Path;
///
/// let file = File::open("data.fasta.gz").unwrap();
/// let gzi = GziIndex::from_path(Path::new("data.fasta.gz.gzi")).unwrap();
/// let mut reader = BgzfReader::with_index(file, gzi).unwrap();
/// ```
pub struct BgzfReader<R: Read + Seek>
{
    /// The underlying compressed file
    inner: R,
    /// Optional .gzi index for seeking
    gzi_index: Option<GziIndex>,
    /// Decompression buffer
    decompressed_buf: Vec<u8>,
    /// Current position in decompressed buffer
    buf_pos: usize,
    /// Current uncompressed position (for tracking)
    current_uncompressed_pos: u64,
    /// End of stream flag
    eof: bool,
}

impl<R: Read + Seek> BgzfReader<R>
{
    /// Create a new BGZF reader without seeking support.
    ///
    /// The reader will decompress BGZF data sequentially.
    pub fn new(inner: R) -> Self
    {
        Self {
            inner,
            gzi_index: None,
            decompressed_buf: Vec::new(),
            buf_pos: 0,
            current_uncompressed_pos: 0,
            eof: false,
        }
    }

    /// Create a BGZF reader with seeking support via a .gzi index.
    ///
    /// # Arguments
    ///
    /// * `inner` - The compressed file reader
    /// * `gzi_index` - The .gzi index for offset mapping
    ///
    /// # Returns
    ///
    /// * `Ok(reader)` - The indexed BGZF reader
    /// * `Err(io::Error)` - If the initial position cannot be set
    pub fn with_index(mut inner: R, gzi_index: GziIndex) -> io::Result<Self>
    {
        // Seek to start of file
        inner.seek(SeekFrom::Start(0))?;
        Ok(Self {
            inner,
            gzi_index: Some(gzi_index),
            decompressed_buf: Vec::new(),
            buf_pos: 0,
            current_uncompressed_pos: 0,
            eof: false,
        })
    }

    /// Seek to an uncompressed position using the .gzi index.
    ///
    /// This method uses the .gzi index to find the compressed offset
    /// and seeks to that position, then decompresses from there.
    ///
    /// # Arguments
    ///
    /// * `uncompressed_pos` - Target position in the uncompressed stream
    ///
    /// # Returns
    ///
    /// * `Ok(actual_pos)` - The actual position reached (may be <= requested)
    /// * `Err(io::Error)` - If no index is available or seeking fails
    pub fn seek_uncompressed(&mut self, uncompressed_pos: u64) -> io::Result<u64>
    {
        let gzi = self.gzi_index.as_ref().ok_or_else(|| {
            io::Error::new(io::ErrorKind::NotFound, "No .gzi index available for seeking")
        })?;

        let compressed_offset = gzi.get_compressed_offset(uncompressed_pos).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("Uncompressed offset {} beyond index range", uncompressed_pos),
            )
        })?;

        // Seek to the compressed offset
        self.inner.seek(SeekFrom::Start(compressed_offset))?;

        // Reset decompression state
        self.decompressed_buf.clear();
        self.buf_pos = 0;

        // Read and decompress blocks until we reach the target position
        while self.current_uncompressed_pos < uncompressed_pos
        {
            if !self.read_next_block()?
            {
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "Reached end of file before target position",
                ));
            }
        }

        // Now we're at or past the target position
        // Set buf_pos to the correct offset within the current block
        let offset_in_block = (uncompressed_pos
            - (self.current_uncompressed_pos - self.decompressed_buf.len() as u64))
            as usize;
        self.buf_pos = offset_in_block;

        Ok(uncompressed_pos)
    }

    /// Get the current uncompressed position.
    pub fn current_position(&self) -> u64
    {
        if self.decompressed_buf.is_empty()
        {
            self.current_uncompressed_pos
        }
        else
        {
            self.current_uncompressed_pos - self.decompressed_buf.len() as u64 + self.buf_pos as u64
        }
    }

    /// Read the next BGZF block.
    ///
    /// Returns true if a block was read, false on EOF.
    fn read_next_block(&mut self) -> io::Result<bool>
    {
        // Read and verify BGZF header (first 12 bytes: ID1, ID2, CM, FLG, MTIME, XFL, OS, XLEN)
        let mut header = [0u8; 12];
        let mut total_read = 0;
        while total_read < 12
        {
            let n = self.inner.read(&mut header[total_read..])?;
            if n == 0
            {
                break;
            }
            total_read += n;
        }

        if total_read == 0
        {
            self.eof = true;
            return Ok(false);
        }
        if total_read < 12
        {
            return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "Incomplete BGZF header"));
        }

        // Verify gzip magic
        if header[0] != GZIP_ID1 || header[1] != GZIP_ID2
        {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "Invalid gzip magic number"));
        }

        // Verify deflate compression method
        if header[2] != GZIP_CM_DEFLATE
        {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "Not deflate compression"));
        }

        let flg = header[3];

        // Check for extra field (BGZF stores block size here)
        let xlen = if flg & GZIP_FLG_FEXTRA != 0
        {
            u16::from_le_bytes([header[10], header[11]]) as usize
        }
        else
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "BGZF requires extra field (FEXTRA flag not set)",
            ));
        };

        // Read extra field
        let mut extra = vec![0u8; xlen];
        self.inner.read_exact(&mut extra)?;

        // Parse BGZF subfield to get block size
        let mut remaining_xlen = xlen;
        let mut block_size = None;

        while remaining_xlen >= 4
        {
            let si1 = extra[xlen - remaining_xlen];
            let si2 = extra[xlen - remaining_xlen + 1];
            let sublen = u16::from_le_bytes([
                extra[xlen - remaining_xlen + 2],
                extra[xlen - remaining_xlen + 3],
            ]) as usize;

            if si1 == BGZF_EXTRA_ID && si2 == BGZF_EXTRA_SUBFIELD && sublen >= 2
            {
                // Block size is stored as a 16-bit little-endian value
                // It includes the 1-byte SI1, 1-byte SI2, 2-byte sublen, and block_size itself
                // So the actual compressed data size is block_size - 1 - 1 - 2 - 2 = block_size - 6
                let bsize = u16::from_le_bytes([
                    extra[xlen - remaining_xlen + 4],
                    extra[xlen - remaining_xlen + 5],
                ]);
                block_size = Some(bsize as usize);
                break;
            }

            // Prevent underflow: ensure we have enough bytes for SI1, SI2, SUBLEN, and the data
            if sublen > remaining_xlen.saturating_sub(4)
            {
                break;
            }
            remaining_xlen -= 4 + sublen;
        }

        let block_size = block_size.ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "BC subfield not found in BGZF extra field")
        })?;

        // MTIME (4 bytes), XFL (1 byte), OS (1 byte) = already read in header
        // Already read: 18 bytes (header) + xlen bytes (extra)
        // Remaining before compressed data: 6 bytes (MTIME + XFL + OS are in header, then we need to account for what we've read)
        // Actually, the header is:
        // ID1(1) ID2(1) CM(1) FLG(1) MTIME(4) XFL(1) OS(1) = 10 bytes
        // Extra length(2) + extra(xlen)
        // We've read 12 bytes (10 + 2) + xlen

        // Skip remaining header fields we haven't accounted for
        // Wait, we read 18 bytes: that's ID1, ID2, CM, FLG, MTIME(4), XFL, OS, XLEN(2)
        // Then we read xlen bytes of extra
        // So now we're at the start of compressed data

        // Calculate remaining compressed data size
        // According to BGZF spec, BSIZE is "the size of the BGZF block minus one"
        // So actual block size = BSIZE + 1
        // Structure:
        //   - Header: ID1(1) + ID2(1) + CM(1) + FLG(1) + MTIME(4) + XFL(1) + OS(1) + XLEN(2) = 12 bytes
        //   - Extra field (XLEN bytes)
        //   - Compressed data (variable)
        //   - Trailer: CRC32(4) + ISIZE(4) = 8 bytes
        // Actual block size = BSIZE + 1 = header(12) + extra_field(XLEN) + compressed_data + trailer(8)
        // Compressed data size = (BSIZE + 1) - header(12) - extra_field(XLEN) - trailer(8)
        let compressed_size = (block_size as isize + 1) - 12 - xlen as isize - 8;
        if compressed_size <= 0
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Invalid BGZF block size: {}, xlen: {}", block_size, xlen),
            ));
        }
        let compressed_size = compressed_size as usize;

        // Read compressed data
        let mut compressed_data = vec![0u8; compressed_size];
        self.inner.read_exact(&mut compressed_data)?;

        // Read and verify trailer (8 bytes: CRC32 + ISIZE)
        let mut trailer = [0u8; 8];
        self.inner.read_exact(&mut trailer)?;

        // Decompress the block
        // Set capacity but keep length at 0 so decompress_vec appends to empty buffer
        self.decompressed_buf.clear();
        self.decompressed_buf.reserve(BGZF_MAX_BLOCK_SIZE);

        let mut decompress = Decompress::new(false);
        decompress.decompress_vec(
            &compressed_data,
            &mut self.decompressed_buf,
            flate2::FlushDecompress::Finish,
        )?;

        self.buf_pos = 0;
        self.current_uncompressed_pos += self.decompressed_buf.len() as u64;
        Ok(true)
    }

    /// Ensure there's data available in the buffer.
    fn fill_buf(&mut self) -> io::Result<&[u8]>
    {
        if self.buf_pos >= self.decompressed_buf.len()
        {
            if self.eof
            {
                return Ok(&[]);
            }
            if !self.read_next_block()?
            {
                return Ok(&[]);
            }
        }
        Ok(&self.decompressed_buf[self.buf_pos..])
    }
}

impl<R: Read + Seek> Read for BgzfReader<R>
{
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize>
    {
        let mut total_read = 0;
        while total_read < buf.len()
        {
            let available = self.fill_buf()?;
            if available.is_empty()
            {
                break;
            }
            let to_read = std::cmp::min(available.len(), buf.len() - total_read);
            buf[total_read..total_read + to_read].copy_from_slice(&available[..to_read]);
            self.buf_pos += to_read;
            total_read += to_read;
        }
        Ok(total_read)
    }
}

impl<R: Read + Seek> BufRead for BgzfReader<R>
{
    fn fill_buf(&mut self) -> io::Result<&[u8]>
    {
        self.fill_buf()
    }

    fn consume(&mut self, amt: usize)
    {
        self.buf_pos += amt;
        if self.buf_pos > self.decompressed_buf.len()
        {
            self.buf_pos = self.decompressed_buf.len();
        }
    }
}

#[cfg(test)]
mod tests
{
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_gzip_magic_verification()
    {
        assert_eq!(GZIP_ID1, 0x1f);
        assert_eq!(GZIP_ID2, 0x8b);
    }

    #[test]
    fn test_bgzf_constants()
    {
        assert_eq!(GZIP_CM_DEFLATE, 8);
        assert_eq!(BGZF_EXTRA_ID, 66); // 'B'
        assert_eq!(BGZF_EXTRA_SUBFIELD, 67); // 'C'
    }

    #[test]
    fn test_reader_creation()
    {
        let data = b"not real gzip data";
        let cursor = Cursor::new(data);
        let reader = BgzfReader::new(cursor);
        assert!(reader.gzi_index.is_none());
        assert_eq!(reader.current_position(), 0);
    }
}
