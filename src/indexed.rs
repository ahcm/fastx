//! Indexed FASTA/FASTQ reader for random access by sequence ID.
//!
//! This module provides `IndexedFastXReader` which enables efficient random access
//! to bgzip-compressed FASTA files using .fai and .gzi indexes.

use crate::bgzf::BgzfReader;
use crate::fai::{FaiEntry, FaiIndex};
use crate::gzi::GziIndex;
use crate::FastX::FastARecord;
use std::fs::File;
use std::io::{self, Read, Seek};
use std::path::Path;

/// An indexed FASTA/FASTQ reader supporting random access by sequence ID.
///
/// This reader uses both .fai (for sequence metadata) and .gzi (for gzip seeking)
/// indexes to efficiently fetch specific sequences without reading the entire file.
///
/// # Type Parameters
///
/// * `R` - The underlying reader type (must implement Read and Seek)
///
/// # Example
///
/// ```no_run
/// use fastx::indexed::IndexedFastXReader;
/// use fastx::FastX::FastXRead;
/// use std::path::Path;
///
/// let mut reader = IndexedFastXReader::from_path(Path::new("data.fasta.gz")).unwrap();
///
/// // Fetch a specific sequence by ID
/// if let Ok(record) = reader.fetch("chr1") {
///     println!("{}: {} bp", record.id(), record.seq_len());
/// }
/// ```
///
/// # URL Support
///
/// With the `url` feature enabled, you can also read from HTTP/HTTPS URLs:
///
/// ```no_run,ignore
/// use fastx::indexed::IndexedFastXReader;
///
/// let mut reader = IndexedFastXReader::from_url(
///     "https://example.com/data.fasta.gz",
///     "https://example.com/data.fasta.gz.fai",
///     "https://example.com/data.fasta.gz.gzi"
/// ).unwrap();
/// ```
pub struct IndexedFastXReader<R: Read + Seek>
{
    /// The BGZF reader for decompression
    reader: BgzfReader<R>,
    /// The FASTA index for sequence lookup
    fai_index: FaiIndex,
}

/// Type alias for local file reading
pub type LocalIndexedFastXReader = IndexedFastXReader<File>;

impl IndexedFastXReader<File>
{
    /// Open an indexed FASTA file from a local path.
    ///
    /// This looks for companion index files (.fai and optionally .gzi) alongside
    /// the specified file.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the FASTA file (can be .fasta.gz or similar)
    ///
    /// # Returns
    ///
    /// * `Ok(reader)` - The indexed reader ready for use
    /// * `Err(io::Error)` - If files cannot be opened or indexes are missing
    ///
    /// # Index Files
    ///
    /// For a file like `data.fasta.gz`:
    /// - `data.fasta.gz.fai` or `data.fasta.fai` - Required FASTA index
    /// - `data.fasta.gz.gzi` or `data.fasta.gzi` - Required gzip index for compressed files
    ///
    /// # Example
    ///
    /// ```no_run
    /// use fastx::indexed::IndexedFastXReader;
    /// use std::path::Path;
    ///
    /// let mut reader = IndexedFastXReader::from_path(Path::new("data.fasta.gz")).unwrap();
    /// ```
    pub fn from_path(path: &Path) -> io::Result<Self>
    {
        // Try to find .fai index
        let fai_path = find_index_file(path, "fai").ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::NotFound,
                format!(
                    "FAI index not found for {} (expected {}.fai or {}.gz.fai)",
                    path.display(),
                    path.with_extension("").display(),
                    path.with_extension("").display(),
                ),
            )
        })?;

        let fai_index = FaiIndex::from_path(&fai_path)?;

        // Check if file is gzip compressed and look for .gzi
        let is_gzip = path.extension().map(|e| e == "gz").unwrap_or(false);

        let file = File::open(path)?;

        let reader = if is_gzip
        {
            // Try to find .gzi index
            if let Some(gzi_path) = find_index_file(path, "gzi")
            {
                let gzi_index = GziIndex::from_path(&gzi_path)?;
                BgzfReader::with_index(file, gzi_index)?
            }
            else
            {
                return Err(io::Error::new(
                    io::ErrorKind::NotFound,
                    format!(
                        "GZI index not found for compressed file {} (expected {}.gzi)",
                        path.display(),
                        path.with_extension("").display()
                    ),
                ));
            }
        }
        else
        {
            return Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "Uncompressed files not yet supported, please use bgzip-compressed files",
            ));
        };

        Ok(Self { reader, fai_index })
    }

    /// Open an indexed FASTA file from HTTP/HTTPS URLs.
    ///
    /// This requires the `url` feature to be enabled.
    ///
    /// # Arguments
    ///
    /// * `data_url` - URL to the FASTA data file (.fasta.gz)
    /// * `fai_url` - URL to the .fai index file
    /// * `gzi_url` - URL to the .gzi index file
    ///
    /// # Returns
    ///
    /// * `Ok(reader)` - The indexed reader ready for use
    /// * `Err(io::Error)` - If URLs are invalid or requests fail
    ///
    /// # Example
    ///
    /// ```no_run,ignore
    /// use fastx::indexed::IndexedFastXReader;
    ///
    /// let mut reader = IndexedFastXReader::from_url(
    ///     "https://example.com/data.fasta.gz",
    ///     "https://example.com/data.fasta.gz.fai",
    ///     "https://example.com/data.fasta.gz.gzi"
    /// ).unwrap();
    ///
    /// let record = reader.fetch("chr1").unwrap();
    /// println!("{}: {} bp", record.id(), record.seq_len());
    /// ```
    #[cfg(feature = "url")]
    pub fn from_url(
        data_url: impl Into<String>,
        fai_url: impl Into<String>,
        gzi_url: impl Into<String>,
    ) -> io::Result<IndexedFastXReader<crate::remote::RemoteReader>>
    {
        use crate::remote::RemoteReader;

        // Fetch and parse the FAI index
        let fai_url = fai_url.into();
        let fai_data = fetch_url(&fai_url)?;
        let fai_index = parse_fai_from_bytes(&fai_data)?;

        // Fetch and parse the GZI index
        let gzi_url = gzi_url.into();
        let gzi_data = fetch_url(&gzi_url)?;
        let gzi_index = parse_gzi_from_bytes(&gzi_data)?;

        // Create the remote reader
        let remote_reader = RemoteReader::new(data_url)?;
        let reader = BgzfReader::with_index(remote_reader, gzi_index)?;

        Ok(IndexedFastXReader { reader, fai_index })
    }
}

impl<R: Read + Seek> IndexedFastXReader<R>
{
    /// Fetch a sequence by its ID.
    ///
    /// Reads the entire sequence from the file using the index.
    ///
    /// # Arguments
    ///
    /// * `seq_id` - The sequence identifier (e.g., "chr1", "gene123")
    ///
    /// # Returns
    ///
    /// * `Ok(FastARecord)` - The fetched sequence record
    /// * `Err(io::Error)` - If the sequence is not found or reading fails
    ///
    /// # Example
    ///
    /// ```no_run
    /// use fastx::indexed::IndexedFastXReader;
    /// use fastx::FastX::FastXRead;
    /// use std::path::Path;
    ///
    /// let mut reader = IndexedFastXReader::from_path(Path::new("data.fasta.gz")).unwrap();
    ///
    /// match reader.fetch("chr1") {
    ///     Ok(record) => println!("Got sequence {}: {} bp", record.id(), record.seq_len()),
    ///     Err(e) => eprintln!("Error: {}", e),
    /// }
    /// ```
    pub fn fetch(&mut self, seq_id: &str) -> io::Result<FastARecord>
    {
        let entry = self.fai_index.get(seq_id).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::NotFound,
                format!("Sequence '{}' not found in index", seq_id),
            )
        })?;

        // Clone the entry to avoid borrowing issues
        let entry = entry.clone();
        self.fetch_entry(&entry)
    }

    /// Fetch a specific region of a sequence.
    ///
    /// # Arguments
    ///
    /// * `seq_id` - The sequence identifier
    /// * `start` - 0-based start position
    /// * `end` - End position (exclusive)
    ///
    /// # Returns
    ///
    /// * `Ok(Vec<u8>)` - The sequence data for the requested region
    /// * `Err(io::Error)` - If the sequence is not found or reading fails
    ///
    /// # Example
    ///
    /// ```no_run
    /// use fastx::indexed::IndexedFastXReader;
    /// use std::path::Path;
    ///
    /// let mut reader = IndexedFastXReader::from_path(Path::new("data.fasta.gz")).unwrap();
    ///
    /// // Fetch bases 1000-2000 of chr1
    /// let region = reader.fetch_range("chr1", 1000, 2000).unwrap();
    /// println!("Region length: {} bp", region.len());
    /// ```
    pub fn fetch_range(&mut self, seq_id: &str, start: u64, end: u64) -> io::Result<Vec<u8>>
    {
        let entry = self.fai_index.get(seq_id).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::NotFound,
                format!("Sequence '{}' not found in index", seq_id),
            )
        })?;

        // Clone to avoid borrowing issues
        let entry = entry.clone();

        if start >= entry.length
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("Start position {} beyond sequence length {}", start, entry.length),
            ));
        }

        let clamped_end = end.min(entry.length);
        let region_length = clamped_end - start;

        // Calculate file offset for start position
        let start_offset = entry.offset_for_position(start);

        // Seek to the start position
        self.reader.seek_uncompressed(start_offset)?;

        // Read the sequence data, handling line wrapping
        let mut seq_data = Vec::with_capacity(region_length as usize);
        let mut remaining = region_length;
        let mut col = start % entry.line_bases;

        while remaining > 0
        {
            // Calculate how much we can read from the current line
            let in_line = std::cmp::min(remaining, entry.line_bases - col);

            // Read that many bytes
            let mut buf = vec![0u8; in_line as usize];
            let n = self.reader.read(&mut buf)?;
            if n == 0
            {
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "Unexpected end of file while reading sequence",
                ));
            }
            seq_data.extend_from_slice(&buf[..n]);
            remaining -= n as u64;

            // Skip the newline
            if col + n as u64 >= entry.line_bases && remaining > 0
            {
                let mut newline = [0u8; 1];
                self.reader.read_exact(&mut newline)?;
                if newline[0] != b'\n'
                {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "Expected newline after sequence line",
                    ));
                }
            }

            col = 0;
        }

        Ok(seq_data)
    }

    /// Fetch a sequence using its FAI entry directly.
    fn fetch_entry(&mut self, entry: &FaiEntry) -> io::Result<FastARecord>
    {
        // Seek to the sequence start
        self.reader.seek_uncompressed(entry.offset)?;

        // Read the header line (first line after seeking)
        let mut header = String::new();
        let mut header_byte = [0u8; 1];
        loop
        {
            let n = self.reader.read(&mut header_byte)?;
            if n == 0 || header_byte[0] == b'\n'
            {
                break;
            }
            header.push(header_byte[0] as char);
        }

        // Read the entire sequence
        let mut raw_seq = Vec::new();
        let mut total_bases = 0u64;
        let mut buf = [0u8; 8192];

        while total_bases < entry.length
        {
            let n = self.reader.read(&mut buf)?;
            if n == 0
            {
                break;
            }

            // Filter out newlines
            for &byte in &buf[..n]
            {
                if byte != b'\n'
                {
                    raw_seq.push(byte);
                    total_bases += 1;
                    if total_bases >= entry.length
                    {
                        break;
                    }
                }
                else
                {
                    // Keep newline in raw_seq for compatibility with FastARecord
                    raw_seq.push(byte);
                }
            }
        }

        Ok(FastARecord {
            name: header,
            raw_seq,
        })
    }

    /// Get a reference to the FAI index.
    pub fn index(&self) -> &FaiIndex
    {
        &self.fai_index
    }

    /// Check if a sequence exists in the index.
    ///
    /// # Arguments
    ///
    /// * `seq_id` - The sequence identifier to check
    pub fn contains(&self, seq_id: &str) -> bool
    {
        self.fai_index.contains(seq_id)
    }

    /// Get all sequence names in the index.
    pub fn sequence_names(&self) -> Vec<&str>
    {
        self.fai_index.sequence_names().collect()
    }
}

/// Fetch data from a URL (requires `url` feature).
#[cfg(feature = "url")]
#[allow(dead_code)]
fn fetch_url(url: &str) -> io::Result<Vec<u8>>
{
    let agent = ureq::Agent::new_with_defaults();

    let response = agent.get(url).call().map_err(|e| {
        io::Error::new(
            io::ErrorKind::ConnectionRefused,
            format!("HTTP GET request failed for {}: {}", url, e),
        )
    })?;

    let data = response.into_body().read_to_vec().map_err(|e| {
        io::Error::new(
            io::ErrorKind::ConnectionRefused,
            format!("Failed to read response body: {}", e),
        )
    })?;

    Ok(data)
}

/// Parse FAI index from bytes (for URL support).
#[allow(dead_code)]
fn parse_fai_from_bytes(data: &[u8]) -> io::Result<FaiIndex>
{
    use crate::fai::FaiEntry;
    use std::collections::HashMap;

    let text = std::str::from_utf8(data)
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "FAI data is not valid UTF-8"))?;

    let mut entries = HashMap::new();

    for (line_num, line) in text.lines().enumerate()
    {
        let line = line.trim();

        // Skip empty lines and comments
        if line.is_empty() || line.starts_with('#')
        {
            continue;
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() != 5
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "Invalid FAI format at line {}: expected 5 fields, got {}",
                    line_num + 1,
                    parts.len()
                ),
            ));
        }

        let name = parts[0].to_string();
        let length = parts[1].parse::<u64>().map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Invalid length at line {}: '{}'", line_num + 1, parts[1]),
            )
        })?;
        let offset = parts[2].parse::<u64>().map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Invalid offset at line {}: '{}'", line_num + 1, parts[2]),
            )
        })?;
        let line_bases = parts[3].parse::<u64>().map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Invalid line_bases at line {}: '{}'", line_num + 1, parts[3]),
            )
        })?;
        let line_width = parts[4].parse::<u64>().map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Invalid line_width at line {}: '{}'", line_num + 1, parts[4]),
            )
        })?;

        if line_width < line_bases
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "Invalid line_width < line_bases at line {}: {} < {}",
                    line_num + 1,
                    line_width,
                    line_bases
                ),
            ));
        }

        let entry = FaiEntry {
            name,
            length,
            offset,
            line_bases,
            line_width,
        };

        entries.insert(entry.name.clone(), entry);
    }

    // Use internal constructor to create FaiIndex
    Ok(FaiIndex { entries })
}

/// Parse GZI index from bytes (for URL support).
#[allow(dead_code)]
fn parse_gzi_from_bytes(data: &[u8]) -> io::Result<GziIndex>
{
    if data.len() < 8
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "GZI data too short (less than 8 bytes)",
        ));
    }

    // Read number of entries (little-endian u64)
    let num_entries = u64::from_le_bytes([
        data[0], data[1], data[2], data[3], data[4], data[5], data[6], data[7],
    ]) as usize;

    let expected_size = 8 + num_entries * 16;
    if data.len() < expected_size
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("GZI data too short: expected {} bytes, got {}", expected_size, data.len()),
        ));
    }

    let mut entries = Vec::with_capacity(num_entries);
    let mut offset = 8;

    for _ in 0..num_entries
    {
        let compressed = u64::from_le_bytes([
            data[offset],
            data[offset + 1],
            data[offset + 2],
            data[offset + 3],
            data[offset + 4],
            data[offset + 5],
            data[offset + 6],
            data[offset + 7],
        ]);
        offset += 8;

        let uncompressed = u64::from_le_bytes([
            data[offset],
            data[offset + 1],
            data[offset + 2],
            data[offset + 3],
            data[offset + 4],
            data[offset + 5],
            data[offset + 6],
            data[offset + 7],
        ]);
        offset += 8;

        entries.push((compressed, uncompressed));
    }

    // Verify entries are sorted by uncompressed offset
    for i in 1..entries.len()
    {
        if entries[i].1 < entries[i - 1].1
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "GZI entries not sorted by uncompressed offset",
            ));
        }
    }

    // Use internal constructor to create GziIndex
    Ok(GziIndex { entries })
}

use std::path::PathBuf;

/// Find an index file for a given data file.
///
/// Tries multiple patterns:
/// - For `data.fasta.gz`: tries `data.fasta.gz.fai` then `data.fasta.fai`
/// - For `data.fasta`: tries `data.fasta.fai`
fn find_index_file(path: &Path, ext: &str) -> Option<PathBuf>
{
    let stem = path.with_extension("");

    // Try path + . + ext (e.g., data.fasta.gz.fai)
    let direct = PathBuf::from(format!("{}.{}", path.display(), ext));
    if direct.exists()
    {
        return Some(direct);
    }

    // Try stem + . + ext (e.g., data.fasta.fai for data.fasta.gz)
    let stem_index = PathBuf::from(format!("{}.{}", stem.display(), ext));
    if stem_index.exists()
    {
        return Some(stem_index);
    }

    None
}

#[cfg(test)]
mod tests
{
    use super::*;

    #[test]
    fn test_find_index_file()
    {
        // Create test files
        let fasta_path = Path::new("test_find.fasta.gz");
        let fai1 = Path::new("test_find.fasta.gz.fai");
        let fai2 = Path::new("test_find.fasta.fai");

        std::fs::write(fasta_path, b">test\nACGT\n").unwrap();

        // Test with .gz.fai extension
        std::fs::write(fai1, b"test\t4\t6\n").unwrap();
        let result = find_index_file(fasta_path, "fai");
        assert!(result.is_some());
        assert_eq!(result.unwrap(), fai1);
        std::fs::remove_file(fai1).unwrap();

        // Test with .fai extension (for .gz file)
        std::fs::write(fai2, b"test\t4\t6\n").unwrap();
        let result = find_index_file(fasta_path, "fai");
        assert!(result.is_some());
        assert_eq!(result.unwrap(), fai2);

        // Cleanup
        std::fs::remove_file(fai2).unwrap();
        std::fs::remove_file(fasta_path).unwrap();
    }

    #[test]
    fn test_index_file_not_found()
    {
        let path = Path::new("nonexistent.fasta.gz");
        let result = find_index_file(path, "fai");
        assert!(result.is_none());
    }

    #[test]
    fn test_parse_fai_from_bytes()
    {
        let data = b"chr1\t1000\t0\t80\t81\nchr2\t2000\t1000\t80\t81\n";
        let index = parse_fai_from_bytes(data).unwrap();
        assert_eq!(index.len(), 2);
        assert!(index.contains("chr1"));
        assert!(index.contains("chr2"));

        let chr1 = index.get("chr1").unwrap();
        assert_eq!(chr1.length, 1000);
        assert_eq!(chr1.offset, 0);
    }

    #[test]
    fn test_parse_gzi_from_bytes()
    {
        let data: Vec<u8> = vec![
            2, 0, 0, 0, 0, 0, 0, 0, // num_entries = 2
            0, 0, 0, 0, 0, 0, 0, 0, // Entry 0: compressed = 0
            0, 0, 0, 0, 0, 0, 0, 0, // Entry 0: uncompressed = 0
            100, 0, 0, 0, 0, 0, 0, 0, // Entry 1: compressed = 100
            0, 100, 0, 0, 0, 0, 0, 0, // Entry 1: uncompressed = 10000
        ];
        let index = parse_gzi_from_bytes(&data).unwrap();
        assert_eq!(index.len(), 2);
        assert_eq!(index.get_compressed_offset(0), Some(0));
        assert_eq!(index.get_compressed_offset(5000), Some(0));
    }
}
