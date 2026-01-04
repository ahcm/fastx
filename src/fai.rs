//! FASTA index (.fai) parser.
//!
//! The .fai format (samtools faidx) stores sequence metadata for random access.
//!
//! Format: tab-separated values per line
//! NAME\tLENGTH\tOFFSET\tLINEBASES\tLINEWIDTH
//!
//! Where:
//! - NAME: Sequence identifier
//! - LENGTH: Total sequence length
//! - OFFSET: Byte offset in uncompressed file where sequence starts
//! - LINEBASES: Number of bases per line
//! - LINEWIDTH: Total bytes per line (including newline)

use std::collections::HashMap;
use std::io;
use std::io::BufRead;
use std::path::Path;

/// A single entry in a FASTA index.
///
/// Contains metadata about a sequence in a FASTA file,
/// enabling efficient random access to that sequence.
///
/// # Fields
///
/// * `name` - Sequence identifier
/// * `length` - Total length of the sequence (bases)
/// * `offset` - Byte offset in the file where this sequence starts
/// * `line_bases` - Number of bases per line in the sequence data
/// * `line_width` - Total bytes per line (bases + newlines)
#[derive(Debug, Clone, PartialEq)]
pub struct FaiEntry
{
    /// Sequence identifier (e.g., "chr1")
    pub name: String,
    /// Total sequence length in bases
    pub length: u64,
    /// Byte offset in uncompressed file where sequence data starts (after header)
    pub offset: u64,
    /// Number of bases per line
    pub line_bases: u64,
    /// Total bytes per line (including newline)
    pub line_width: u64,
}

impl FaiEntry
{
    /// Calculate the byte offset for a specific region within this sequence.
    ///
    /// Given a start position (0-based), calculates where in the file to start reading.
    /// Accounts for line wrapping in the FASTA file.
    ///
    /// # Arguments
    ///
    /// * `start` - 0-based start position within the sequence
    ///
    /// # Returns
    ///
    /// The byte offset from the beginning of the file
    ///
    /// # Example
    ///
    /// ```
    /// use fastx::fai::FaiEntry;
    ///
    /// // Sequence with 80 bases per line, 81 bytes per line (including newline)
    /// let entry = FaiEntry {
    ///     name: "chr1".to_string(),
    ///     length: 1000,
    ///     offset: 100,
    ///     line_bases: 80,
    ///     line_width: 81,
    /// };
    ///
    /// // Position 100 is on line 2 (0-based), column 20
    /// // Line 0: offsets 0-79
    /// // Line 1: offsets 80-159 (position 100 is here, at offset 20 from line start)
    /// // File offset = sequence_offset + (1 * line_width) + 20
    /// let file_offset = entry.offset_for_position(100);
    /// assert_eq!(file_offset, 100 + 81 + 20);
    /// ```
    pub fn offset_for_position(&self, start: u64) -> u64
    {
        let full_lines = start / self.line_bases;
        let col = start % self.line_bases;
        self.offset + (full_lines * self.line_width) + col
    }

    /// Calculate the length of a region, accounting for line wrapping.
    ///
    /// Returns the number of sequence bases in the specified region,
    /// clamped to the sequence length.
    ///
    /// # Arguments
    ///
    /// * `start` - 0-based start position
    /// * `end` - End position (exclusive)
    pub fn region_length(&self, start: u64, end: u64) -> u64
    {
        let clamped_end = end.min(self.length);
        let clamped_start = start.min(self.length);
        clamped_end.saturating_sub(clamped_start)
    }
}

/// A FASTA index for random access to sequences.
///
/// The .fai format is used by samtools faidx and other tools to enable
/// efficient random access to FASTA files.
///
/// # Example
///
/// ```no_run
/// use fastx::fai::FaiIndex;
/// use std::path::Path;
///
/// let index = FaiIndex::from_path(Path::new("data.fasta.fai")).unwrap();
/// if let Some(entry) = index.get("chr1") {
///     println!("chr1 length: {}", entry.length);
/// }
/// ```
#[derive(Debug, Clone)]
pub struct FaiIndex
{
    pub entries: HashMap<String, FaiEntry>,
}

impl FaiIndex
{
    /// Load a .fai index from a file.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the .fai index file
    ///
    /// # Returns
    ///
    /// * `Ok(FaiIndex)` - The parsed index
    /// * `Err(io::Error)` - If the file cannot be read or the format is invalid
    ///
    /// # Example
    ///
    /// ```no_run
    /// use fastx::fai::FaiIndex;
    /// use std::path::Path;
    ///
    /// let index = FaiIndex::from_path(Path::new("sequences.fasta.fai")).unwrap();
    /// ```
    pub fn from_path(path: &Path) -> io::Result<Self>
    {
        let file = std::fs::File::open(path)?;
        let reader = io::BufReader::new(file);
        let mut entries = HashMap::new();

        for (line_num, line_result) in reader.lines().enumerate()
        {
            let line = line_result?;
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

            // Use the name as key (first column)
            entries.insert(entry.name.clone(), entry);
        }

        Ok(FaiIndex { entries })
    }

    /// Get an entry by sequence name.
    ///
    /// # Arguments
    ///
    /// * `name` - The sequence identifier
    ///
    /// # Returns
    ///
    /// * `Some(&FaiEntry)` - If the sequence exists in the index
    /// * `None` - If the sequence is not found
    ///
    /// # Example
    ///
    /// ```no_run
    /// use fastx::fai::FaiIndex;
    /// use std::path::Path;
    ///
    /// let index = FaiIndex::from_path(Path::new("data.fasta.fai")).unwrap();
    /// if let Some(entry) = index.get("chr1") {
    ///     println!("chr1 length: {}", entry.length);
    /// }
    /// ```
    pub fn get(&self, name: &str) -> Option<&FaiEntry>
    {
        self.entries.get(name)
    }

    /// Check if a sequence exists in the index.
    ///
    /// # Arguments
    ///
    /// * `name` - The sequence identifier
    pub fn contains(&self, name: &str) -> bool
    {
        self.entries.contains_key(name)
    }

    /// Get the number of sequences in the index.
    pub fn len(&self) -> usize
    {
        self.entries.len()
    }

    /// Check if the index is empty.
    pub fn is_empty(&self) -> bool
    {
        self.entries.is_empty()
    }

    /// Get an iterator over all sequence names in the index.
    pub fn sequence_names(&self) -> impl Iterator<Item = &str>
    {
        self.entries.keys().map(|s| s.as_str())
    }

    /// Get all entries in the index.
    pub fn entries(&self) -> impl Iterator<Item = &FaiEntry>
    {
        self.entries.values()
    }
}

#[cfg(test)]
mod tests
{
    use super::*;

    const TEST_FAI: &str = "chr1\t248956422\t6\t80\t81
chr2\t242193529\t250000000\t80\t81
chr3\t198295559\t493000000\t80\t81
";

    #[test]
    fn test_fai_parsing()
    {
        let path = Path::new("test.fasta.fai");
        std::fs::write(path, TEST_FAI).unwrap();

        let index = FaiIndex::from_path(path).unwrap();

        assert_eq!(index.len(), 3);
        assert!(index.contains("chr1"));
        assert!(index.contains("chr2"));
        assert!(index.contains("chr3"));
        assert!(!index.contains("chr4"));

        let chr1 = index.get("chr1").unwrap();
        assert_eq!(chr1.name, "chr1");
        assert_eq!(chr1.length, 248956422);
        assert_eq!(chr1.offset, 6);
        assert_eq!(chr1.line_bases, 80);
        assert_eq!(chr1.line_width, 81);

        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_offset_for_position()
    {
        let entry = FaiEntry {
            name: "test".to_string(),
            length: 1000,
            offset: 100,
            line_bases: 80,
            line_width: 81,
        };

        // Position 0 -> offset 100
        assert_eq!(entry.offset_for_position(0), 100);

        // Position 79 -> end of first line
        assert_eq!(entry.offset_for_position(79), 179);

        // Position 80 -> start of second line (line 0: 100-180, line 1: 181-...)
        assert_eq!(entry.offset_for_position(80), 181);

        // Position 100 -> line 1, column 20
        // offset = 100 + (1 * 81) + 20 = 201
        assert_eq!(entry.offset_for_position(100), 201);

        // Position 160 -> line 2, column 0
        // offset = 100 + (2 * 81) + 0 = 262
        assert_eq!(entry.offset_for_position(160), 262);
    }

    #[test]
    fn test_region_length()
    {
        let entry = FaiEntry {
            name: "test".to_string(),
            length: 1000,
            offset: 0,
            line_bases: 80,
            line_width: 81,
        };

        // Normal range
        assert_eq!(entry.region_length(100, 200), 100);

        // Clamped to sequence length
        assert_eq!(entry.region_length(900, 2000), 100);

        // Empty range
        assert_eq!(entry.region_length(100, 100), 0);

        // Start beyond sequence
        assert_eq!(entry.region_length(2000, 2000), 0);
    }

    #[test]
    fn test_empty_lines_and_comments()
    {
        let data = "# This is a comment\n\nchr1\t100\t0\t80\t81\n\n";
        let path = Path::new("test_comment.fai");
        std::fs::write(path, data).unwrap();

        let index = FaiIndex::from_path(path).unwrap();
        assert_eq!(index.len(), 1);
        assert!(index.contains("chr1"));

        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_invalid_format()
    {
        // Not enough fields
        let data = "chr1\t100\t0\t80\n";
        let path = Path::new("test_invalid.fai");
        std::fs::write(path, data).unwrap();

        assert!(FaiIndex::from_path(path).is_err());

        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_invalid_line_width()
    {
        // line_width < line_bases
        let data = "chr1\t100\t0\t80\t70\n";
        let path = Path::new("test_invalid_width.fai");
        std::fs::write(path, data).unwrap();

        assert!(FaiIndex::from_path(path).is_err());

        std::fs::remove_file(path).unwrap();
    }
}
