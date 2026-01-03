//! Indexed FASTA/FASTQ reader for random access by sequence ID.
//!
//! This module provides `IndexedFastXReader` which enables efficient random access
//! to bgzip-compressed FASTA files using .fai and .gzi indexes.

use crate::bgzf::BgzfReader;
use crate::fai::{FaiEntry, FaiIndex};
use crate::gzi::GziIndex;
use crate::FastX::FastARecord;
use std::fs::File;
use std::io::{self, Read};
use std::path::Path;

/// An indexed FASTA/FASTQ reader supporting random access by sequence ID.
///
/// This reader uses both .fai (for sequence metadata) and .gzi (for gzip seeking)
/// indexes to efficiently fetch specific sequences without reading the entire file.
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
pub struct IndexedFastXReader
{
    /// The BGZF reader for decompression
    reader: BgzfReader<File>,
    /// The FASTA index for sequence lookup
    fai_index: FaiIndex,
}

impl IndexedFastXReader
{
    /// Open an indexed FASTA file.
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
    /// - `data.fasta.gz.gzi` or `data.fasta.gzi` - Optional gzip index for compressed files
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
        let is_gzip = path
            .extension()
            .map(|e| e == "gz")
            .unwrap_or(false);

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
            // For uncompressed files, we could use a simpler reader
            // For now, require BGZF format
            return Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "Uncompressed files not yet supported, please use bgzip-compressed files",
            ));
        };

        Ok(Self {
            reader,
            fai_index,
        })
    }

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
                format!(
                    "Start position {} beyond sequence length {}",
                    start, entry.length
                ),
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

/// Find an index file for a given data file.
///
/// Tries multiple patterns:
/// - For `data.fasta.gz`: tries `data.fasta.gz.fai` then `data.fasta.fai`
/// - For `data.fasta`: tries `data.fasta.fai`
fn find_index_file(path: &Path, ext: &str) -> Option<PathBuf>
{
    use std::path::PathBuf;

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

use std::path::PathBuf;

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
}
