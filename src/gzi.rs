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
    pub entries: Vec<(u64, u64)>,
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
    pub fn from_path(path: &Path) -> io::Result<Self>
    {
        let mut file = std::fs::File::open(path)?;
        let mut buffer = Vec::new();
        file.read_to_end(&mut buffer)?;
        Self::from_bytes(&buffer)
    }

    /// Load a .gzi index from raw bytes.
    ///
    /// # Arguments
    ///
    /// * `buffer` - Byte slice containing .gzi formatted data
    ///
    /// # Returns
    ///
    /// * `Ok(GziIndex)` - The parsed index
    /// * `Err(io::Error)` - If the data is too short or the format is invalid
    pub fn from_bytes(buffer: &[u8]) -> io::Result<Self>
    {
        if buffer.len() < 8
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "GZI data too short (less than 8 bytes)",
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
                    "GZI data too short: expected {} bytes, got {}",
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
    pub fn get_compressed_offset(&self, uncompressed_offset: u64) -> Option<u64>
    {
        if self.entries.is_empty()
        {
            return None;
        }

        let result = self
            .entries
            .binary_search_by(|(_, unc)| unc.cmp(&uncompressed_offset));

        match result
        {
            Ok(i) => Some(self.entries[i].0),
            Err(0) => Some(self.entries[0].0),
            Err(i) if i >= self.entries.len() => Some(self.entries.last()?.0),
            Err(i) => Some(self.entries[i - 1].0),
        }
    }

    /// Get the uncompressed offset for a given compressed position.
    pub fn get_uncompressed_offset(&self, compressed_offset: u64) -> Option<u64>
    {
        if self.entries.is_empty()
        {
            return None;
        }

        let result = self
            .entries
            .binary_search_by(|(comp, _)| comp.cmp(&compressed_offset));

        match result
        {
            Ok(i) => Some(self.entries[i].1),
            Err(0) => None,
            Err(i) if i >= self.entries.len() => Some(self.entries.last()?.1),
            Err(i) => Some(self.entries[i - 1].1),
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

    #[test]
    fn test_from_bytes()
    {
        let data: Vec<u8> = vec![
            1, 0, 0, 0, 0, 0, 0, 0, // 1 entry
            100, 0, 0, 0, 0, 0, 0, 0, // comp = 100
            0, 0, 0, 0, 0, 0, 0, 0, // uncomp = 0
        ];
        let index = GziIndex::from_bytes(&data).unwrap();
        assert_eq!(index.len(), 1);
        assert_eq!(index.get_compressed_offset(0), Some(100));
    }
}
