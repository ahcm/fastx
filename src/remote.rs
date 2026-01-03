//! Remote file reader with HTTP range request support and block-level caching.
//!
//! This module provides `RemoteReader` which implements `Read` and `Seek`
//! for HTTP/HTTPS URLs using range requests, with intelligent caching to
//! minimize network requests.

#![cfg(feature = "url")]

use std::collections::HashMap;
use std::io::{self, Read, Seek, SeekFrom};
use std::sync::Mutex;
use ureq::Agent;

/// Default block size for caching (64KB).
const DEFAULT_BLOCK_SIZE: u64 = 64 * 1024;

/// A cached block of data from the remote file.
#[derive(Debug, Clone)]
struct CachedBlock
{
    /// Starting offset of this block in the file
    #[allow(dead_code)]
    offset: u64,
    /// The cached data
    data: Vec<u8>,
}

/// A remote file reader with HTTP range request support and caching.
///
/// This reader fetches data from HTTP/HTTPS URLs on demand, caching blocks
/// to minimize network traffic. It implements `Read` and `Seek` for random
/// access to remote files.
///
/// # Caching
///
/// The reader caches 64KB blocks. When data is requested, it fetches the
/// entire block containing that position, serving subsequent reads from
/// the same range from the cache.
///
/// # Example
///
/// ```no_run
/// use fastx::remote::RemoteReader;
///
/// let reader = RemoteReader::new("https://example.com/data.fasta.gz").unwrap();
/// ```
pub struct RemoteReader
{
    /// The base URL
    url: String,
    /// The HTTP agent for making requests
    agent: Agent,
    /// Cache of fetched blocks (offset -> data)
    cache: Mutex<HashMap<u64, CachedBlock>>,
    /// Current position in the file
    pos: u64,
    /// Total file size (cached after first request)
    file_size: Option<u64>,
    /// Block size for caching
    block_size: u64,
}

impl RemoteReader
{
    /// Create a new remote reader for the given URL.
    ///
    /// # Arguments
    ///
    /// * `url` - The HTTP/HTTPS URL to read from
    ///
    /// # Returns
    ///
    /// * `Ok(reader)` - The remote reader ready for use
    /// * `Err(io::Error)` - If the URL is invalid or the request fails
    ///
    /// # Example
    ///
    /// ```no_run
    /// use fastx::remote::RemoteReader;
    ///
    /// let reader = RemoteReader::new("https://example.com/data.fasta.gz").unwrap();
    /// ```
    pub fn new(url: impl Into<String>) -> io::Result<Self>
    {
        let url = url.into();
        let agent = Agent::new_with_defaults();

        // Probe for file size using a HEAD request
        let file_size = Self::get_file_size_for_url(&agent, &url)?;

        Ok(Self {
            url,
            agent,
            cache: Mutex::new(HashMap::new()),
            pos: 0,
            file_size: Some(file_size),
            block_size: DEFAULT_BLOCK_SIZE,
        })
    }

    /// Get the total file size for a URL (static helper).
    fn get_file_size_for_url(agent: &Agent, url: &str) -> io::Result<u64>
    {
        let response = agent.head(url).call().map_err(|e| {
            io::Error::new(
                io::ErrorKind::ConnectionRefused,
                format!("HTTP HEAD request failed: {}", e),
            )
        })?;

        let content_length = response
            .headers()
            .get("Content-Length")
            .and_then(|v| v.to_str().ok())
            .and_then(|s| s.parse::<u64>().ok())
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Missing or invalid Content-Length header",
                )
            })?;

        Ok(content_length)
    }

    /// Set the block size for caching.
    ///
    /// Larger blocks reduce the number of HTTP requests but use more memory.
    ///
    /// # Arguments
    ///
    /// * `size` - Block size in bytes
    pub fn with_block_size(mut self, size: u64) -> Self
    {
        self.block_size = size;
        self
    }

    /// Get the total file size.
    ///
    /// Makes a HEAD request to determine Content-Length if not already cached.
    fn get_file_size(&self) -> io::Result<u64>
    {
        if let Some(size) = self.file_size
        {
            return Ok(size);
        }

        let response = self.agent.head(&self.url).call().map_err(|e| {
            io::Error::new(
                io::ErrorKind::ConnectionRefused,
                format!("HTTP HEAD request failed: {}", e),
            )
        })?;

        let content_length = response
            .headers()
            .get("Content-Length")
            .and_then(|v| v.to_str().ok())
            .and_then(|s| s.parse::<u64>().ok())
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Missing or invalid Content-Length header",
                )
            })?;

        Ok(content_length)
    }

    /// Get the starting offset of the block containing a given position.
    fn block_start(&self, pos: u64) -> u64
    {
        (pos / self.block_size) * self.block_size
    }

    /// Fetch a block from the remote server.
    ///
    /// # Arguments
    ///
    /// * `offset` - Starting offset of the block
    fn fetch_block(&self, offset: u64) -> io::Result<CachedBlock>
    {
        let file_size = self.get_file_size()?;
        let end = std::cmp::min(offset + self.block_size - 1, file_size.saturating_sub(1));

        let range = if offset >= file_size
        {
            return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "Seek beyond end of file"));
        }
        else if end < offset
        {
            // Empty file or offset at end
            format!("bytes={0}-", offset)
        }
        else
        {
            format!("bytes={}-{}", offset, end)
        };

        let response = self
            .agent
            .get(&self.url)
            .header("Range", &range)
            .call()
            .map_err(|e| {
                io::Error::new(
                    io::ErrorKind::ConnectionRefused,
                    format!("HTTP GET request failed: {}", e),
                )
            })?;

        // Check for partial content or OK status
        let status = response.status();
        if status != 206 && status != 200
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Unexpected HTTP status: {}", status),
            ));
        }

        let data = response.into_body().read_to_vec().map_err(|e| {
            io::Error::new(
                io::ErrorKind::ConnectionRefused,
                format!("Failed to read response body: {}", e),
            )
        })?;

        Ok(CachedBlock { offset, data })
    }

    /// Get data at a specific offset, using cache if available.
    ///
    /// # Arguments
    ///
    /// * `offset` - Position in the file to read from
    ///
    /// # Returns
    ///
    /// A slice containing the cached block data
    fn get_data_at(&self, offset: u64) -> io::Result<Vec<u8>>
    {
        let block_start = self.block_start(offset);

        // Check if we need to fetch the block
        if !self
            .cache
            .lock()
            .map_err(|_| io::Error::new(io::ErrorKind::Other, "Cache lock poisoned"))?
            .contains_key(&block_start)
        {
            // Fetch the block
            let block = self.fetch_block(block_start)?;
            let mut cache = self
                .cache
                .lock()
                .map_err(|_| io::Error::new(io::ErrorKind::Other, "Cache lock poisoned"))?;
            cache.insert(block_start, block);
        }

        // Get the data from cache
        let cache = self
            .cache
            .lock()
            .map_err(|_| io::Error::new(io::ErrorKind::Other, "Cache lock poisoned"))?;
        let block = cache.get(&block_start).unwrap();
        let offset_in_block = (offset - block_start) as usize;
        Ok(block.data[offset_in_block..].to_vec())
    }
}

impl Read for RemoteReader
{
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize>
    {
        let file_size = self.get_file_size()?;
        if self.pos >= file_size
        {
            return Ok(0);
        }

        let remaining = file_size - self.pos;
        let to_read = std::cmp::min(buf.len() as u64, remaining) as usize;

        // Fetch data for current position
        let data = self.get_data_at(self.pos)?;

        let actual_read = std::cmp::min(to_read, data.len());
        buf[..actual_read].copy_from_slice(&data[..actual_read]);
        self.pos += actual_read as u64;

        Ok(actual_read)
    }
}

impl Seek for RemoteReader
{
    fn seek(&mut self, pos: SeekFrom) -> io::Result<u64>
    {
        let file_size = self.get_file_size().ok();

        self.pos = match pos
        {
            SeekFrom::Start(n) => n,
            SeekFrom::End(offset) =>
            {
                let size = file_size
                    .ok_or_else(|| io::Error::new(io::ErrorKind::NotFound, "Unknown file size"))?;
                let offset_i64 = offset as i64;
                if offset_i64 < 0
                {
                    size.checked_sub(offset_i64.unsigned_abs()).ok_or_else(|| {
                        io::Error::new(io::ErrorKind::InvalidInput, "Seek before file start")
                    })?
                }
                else
                {
                    size.checked_add(offset as u64).ok_or_else(|| {
                        io::Error::new(io::ErrorKind::InvalidInput, "Seek overflow")
                    })?
                }
            }
            SeekFrom::Current(offset) =>
            {
                let offset_i64 = offset as i64;
                if offset_i64 < 0
                {
                    self.pos
                        .checked_sub(offset_i64.unsigned_abs())
                        .ok_or_else(|| {
                            io::Error::new(io::ErrorKind::InvalidInput, "Seek before file start")
                        })?
                }
                else
                {
                    self.pos.checked_add(offset as u64).ok_or_else(|| {
                        io::Error::new(io::ErrorKind::InvalidInput, "Seek overflow")
                    })?
                }
            }
        };

        Ok(self.pos)
    }
}

#[cfg(test)]
mod tests
{
    use super::*;

    #[test]
    fn test_block_start()
    {
        // Test block_start calculation without making HTTP requests
        let url = "http://example.com/test";
        let agent = Agent::new_with_defaults();

        // Create a reader without probing file size
        let reader = RemoteReader {
            url: url.to_string(),
            agent,
            cache: Mutex::new(HashMap::new()),
            pos: 0,
            file_size: None,
            block_size: DEFAULT_BLOCK_SIZE,
        };

        assert_eq!(reader.block_start(0), 0);
        assert_eq!(reader.block_start(100), 0);
        assert_eq!(reader.block_start(65536), 65536);
        assert_eq!(reader.block_start(70000), 65536);
        assert_eq!(reader.block_start(131072), 131072);
    }

    #[test]
    fn test_block_start_custom_size()
    {
        let url = "http://example.com/test";
        let agent = Agent::new_with_defaults();

        let reader = RemoteReader {
            url: url.to_string(),
            agent,
            cache: Mutex::new(HashMap::new()),
            pos: 0,
            file_size: None,
            block_size: 1024,
        };

        assert_eq!(reader.block_start(0), 0);
        assert_eq!(reader.block_start(100), 0);
        assert_eq!(reader.block_start(1024), 1024);
        assert_eq!(reader.block_start(2000), 1024);
    }
}
