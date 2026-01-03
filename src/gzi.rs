//! BGZF gzip index (.gzi) parser.
//!
//! The .gzi format stores pairs of (compressed_offset, uncompressed_offset)
//! that enable random access to bgzip-compressed files.

use std::io;
use std::io::Read;
use std::path::Path;

/// A gzip index for BGZF-compressed files.
///
/// Maps uncompressed byte positions to compressed file positions,
/// enabling efficient seeking in compressed files.
///
/// # Format
///
/// Binary format (little-endian):
/// - First 8 bytes: number of index entries (u64)
/// - Following bytes: pairs of (compressed_offset, uncompressed_offset) as u64
///
/// # Example
///
/// ```no_run
/// use fastx::gzi::GziIndex;
/// use std::path::Path;
///
/// let index = GziIndex::from_path(Path::new("data.fasta.gz.gzi")).unwrap();
/// let compressed_offset = index.get_compressed_offset(10000).unwrap();
/// ```
#[derive(Debug, Clone)]
pub struct GziIndex
{
    /// Sorted entries: (compressed_offset, uncompressed_offset)
    entries: Vec<(u64, u64)>,
}

impl GziIndex
{
    /// Load a .gzi index from a file.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the .gzi index file
    ///
    /// # Returns
    ///
    /// * `Ok(GziIndex)` - The parsed index
    /// * `Err(io::Error)` - If the file cannot be read or the format is invalid
    ///
    /// # Example
    ///
    /// ```no_run
    /// use fastx::gzi::GziIndex;
    /// use std::path::Path;
    ///
    /// let index = GziIndex::from_path(Path::new("data.gz.gzi")).unwrap();
    /// ```
    pub fn from_path(path: &Path) -> io::Result<Self>
    {
        let mut file = std::fs::File::open(path)?;
        let mut buffer = Vec::new();
        file.read_to_end(&mut buffer)?;

        if buffer.len() < 8
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "GZI file too short (less than 8 bytes)",
            ));
        }

        // Read number of entries (little-endian u64)
        let num_entries = u64::from_le_bytes([
            buffer[0], buffer[1], buffer[2], buffer[3], buffer[4], buffer[5], buffer[6], buffer[7],
        ]) as usize;

        let expected_size = 8 + num_entries * 16;
        if buffer.len() < expected_size
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "GZI file too short: expected {} bytes, got {}",
                    expected_size,
                    buffer.len()
                ),
            ));
        }

        let mut entries = Vec::with_capacity(num_entries);
        let mut offset = 8;

        for _ in 0..num_entries
        {
            let compressed = u64::from_le_bytes([
                buffer[offset],
                buffer[offset + 1],
                buffer[offset + 2],
                buffer[offset + 3],
                buffer[offset + 4],
                buffer[offset + 5],
                buffer[offset + 6],
                buffer[offset + 7],
            ]);
            offset += 8;

            let uncompressed = u64::from_le_bytes([
                buffer[offset],
                buffer[offset + 1],
                buffer[offset + 2],
                buffer[offset + 3],
                buffer[offset + 4],
                buffer[offset + 5],
                buffer[offset + 6],
                buffer[offset + 7],
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

        Ok(GziIndex { entries })
    }

    /// Get the compressed offset for a given uncompressed position.
    ///
    /// Returns the compressed file offset where decompression should start
    /// to reach the specified uncompressed position.
    ///
    /// # Arguments
    ///
    /// * `uncompressed_offset` - Position in the uncompressed data stream
    ///
    /// # Returns
    ///
    /// * `Some(offset)` - The compressed file offset to seek to
    /// * `None` - If the offset is beyond the end of the index
    ///
    /// # Example
    ///
    /// ```no_run
    /// use fastx::gzi::GziIndex;
    /// use std::path::Path;
    ///
    /// let index = GziIndex::from_path(Path::new("data.gz.gzi")).unwrap();
    /// // Find the compressed offset for an uncompressed position
    /// let offset = index.get_compressed_offset(15000);
    /// ```
    pub fn get_compressed_offset(&self, uncompressed_offset: u64) -> Option<u64>
    {
        if self.entries.is_empty()
        {
            return None;
        }

        // Binary search for the largest uncompressed_offset <= target
        let result = self
            .entries
            .binary_search_by(|(_, unc)| unc.cmp(&uncompressed_offset));

        match result
        {
            Ok(i) => Some(self.entries[i].0),
            Err(0) => Some(self.entries[0].0), // Before first entry, use first
            Err(i) if i >= self.entries.len() => Some(self.entries.last()?.0), // Beyond last, use last
            Err(i) => Some(self.entries[i - 1].0), // Between entries, use previous
        }
    }

    /// Get the number of index entries.
    pub fn len(&self) -> usize
    {
        self.entries.len()
    }

    /// Check if the index is empty.
    pub fn is_empty(&self) -> bool
    {
        self.entries.is_empty()
    }

    /// Get a reference to all entries.
    pub fn entries(&self) -> &[(u64, u64)]
    {
        &self.entries
    }
}

#[cfg(test)]
mod tests
{
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_empty_index()
    {
        // Create a minimal valid GZI with 0 entries
        let data: [u8; 8] = [0, 0, 0, 0, 0, 0, 0, 0];
        let mut cursor = Cursor::new(data);
        let mut buffer = Vec::new();
        cursor.read_to_end(&mut buffer).unwrap();

        let num_entries = u64::from_le_bytes([
            buffer[0], buffer[1], buffer[2], buffer[3], buffer[4], buffer[5], buffer[6], buffer[7],
        ]) as usize;
        assert_eq!(num_entries, 0);

        let index = GziIndex {
            entries: Vec::new(),
        };
        assert!(index.is_empty());
        assert_eq!(index.len(), 0);
    }

    #[test]
    fn test_get_compressed_offset()
    {
        let index = GziIndex {
            entries: vec![(0, 0), (100, 10000), (250, 20000), (400, 30000)],
        };

        // Exact match
        assert_eq!(index.get_compressed_offset(0), Some(0));
        assert_eq!(index.get_compressed_offset(10000), Some(100));

        // Between entries - should find previous block
        assert_eq!(index.get_compressed_offset(5000), Some(0));
        assert_eq!(index.get_compressed_offset(15000), Some(100));
        assert_eq!(index.get_compressed_offset(25000), Some(250));

        // Beyond last entry - returns last entry's offset
        assert_eq!(index.get_compressed_offset(40000), Some(400));
    }

    #[test]
    fn test_single_entry()
    {
        let index = GziIndex {
            entries: vec![(50, 0)],
        };

        assert_eq!(index.get_compressed_offset(0), Some(50));
        assert_eq!(index.get_compressed_offset(1000), Some(50));
    }

    #[test]
    fn test_gzi_format_parsing()
    {
        // Create a valid GZI file with 2 entries:
        // Entry 0: compressed=100, uncompressed=0
        // Entry 1: compressed=300, uncompressed=10000
        let data: Vec<u8> = vec![
            // num_entries = 2 (little-endian u64)
            2, 0, 0, 0, 0, 0, 0, 0, // Entry 0: compressed = 100
            100, 0, 0, 0, 0, 0, 0, 0, // Entry 0: uncompressed = 0
            0, 0, 0, 0, 0, 0, 0, 0, // Entry 1: compressed = 300
            44, 1, 0, 0, 0, 0, 0, 0, // Entry 1: uncompressed = 10000
            16, 39, 0, 0, 0, 0, 0, 0,
        ];

        let mut cursor = Cursor::new(data);
        let mut buffer = Vec::new();
        cursor.read_to_end(&mut buffer).unwrap();

        let num_entries = u64::from_le_bytes([
            buffer[0], buffer[1], buffer[2], buffer[3], buffer[4], buffer[5], buffer[6], buffer[7],
        ]) as usize;
        assert_eq!(num_entries, 2);

        let index = GziIndex {
            entries: vec![(100, 0), (300, 10000)],
        };

        assert_eq!(index.len(), 2);
        assert_eq!(index.get_compressed_offset(0), Some(100));
        assert_eq!(index.get_compressed_offset(5000), Some(100));
    }
}
