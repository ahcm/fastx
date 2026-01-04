/*
inspired by seq_io, https://github.com/markschl/seq_io
Copyright (c) 2021 Andreas Hauser <Andreas.Hauser@LMU.de>
License: Attribution-ShareAlike 4.0 International
 */

//! Low-overhead readers for FASTA and FASTQ sequence files.
//!
//! FastX provides efficient parsing of FASTA and FASTQ formatted files, which are
//! commonly used in bioinformatics to store nucleotide or protein sequences.
//!
//! # Features
//!
//! - Zero-copy parsing where possible
//! - Support for gzip-compressed files (.gz)
//! - Iterator-based API for easy processing
//! - Manual read API for fine-grained control
//! - Automatic format detection
//!
//! # Format Overview
//!
//! **FASTA format**: Used to represent nucleotide or peptide sequences.
//! Each record starts with `>` followed by a name/description line, then sequence data.
//!
//! ```text
//! >sequence_name description
//! AGCTTAGCTAGCTACGATCG
//! ```
//!
//! **FASTQ format**: Like FASTA but includes quality scores for each base.
//! Each record has four lines: name (starting with `@`), sequence, separator (`+`), and quality.
//!
//! ```text
//! @sequence_name
//! AGCTTAGCTAGCTACGATCG
//! +
//! !''*((((***+))%%%++
//! ```
//!
//! # Examples
//!
//! ## High-performance reading (recommended for speed)
//!
//! ```no_run
//! use std::io::BufReader;
//! use std::fs::File;
//! use fastx::FastX::{fasta_for_each, FastXRead};
//!
//! let reader = BufReader::new(File::open("sequences.fasta").unwrap());
//! fasta_for_each(reader, |record| {
//!     println!("{} - length: {}", record.id(), record.seq_len());
//! }).unwrap();
//! ```
//!
//! ## Iterator-based reading (convenient)
//!
//! ```no_run
//! use std::io::BufReader;
//! use std::fs::File;
//! use fastx::FastX::{fasta_iter, FastXRead};
//!
//! let reader = BufReader::new(File::open("sequences.fasta").unwrap());
//! for result in fasta_iter(reader) {
//!     let record = result.unwrap();
//!     println!("{} - length: {}", record.id(), record.seq_len());
//! }
//! ```
//!
//! ## Manual reading
//!
//! ```no_run
//! use std::io::BufRead;
//! use fastx::FastX::{FastARecord, FastXRead, reader_from_path};
//! use std::path::Path;
//!
//! let mut reader = reader_from_path(Path::new("sequences.fasta")).unwrap();
//! let mut record = FastARecord::default();
//!
//! while record.read(&mut *reader).unwrap() > 0 {
//!     println!("{}", record);
//! }
//! ```
//!
//! ## Working with gzip-compressed files
//!
//! ```no_run
//! use fastx::FastX::{reader_from_path, fastq_iter, FastXRead};
//! use std::path::Path;
//!
//! // Automatically detects .gz extension and decompresses
//! let reader = reader_from_path(Path::new("sequences.fastq.gz")).unwrap();
//! for result in fastq_iter(reader) {
//!     let record = result.unwrap();
//!     println!("{}", record.id());
//! }
//! ```
//!
//! ## Random access with indexed files
//!
//! ```no_run
//! use fastx::indexed::IndexedFastXReader;
//! use fastx::FastX::FastXRead;
//! use std::path::Path;
//!
//! // Open an indexed FASTA file (requires .fai and .gzi files)
//! let mut reader = IndexedFastXReader::from_path(Path::new("data.fasta.gz")).unwrap();
//!
//! // Fetch a specific sequence by ID
//! let record = reader.fetch("chr1").unwrap();
//! println!("{}: {} bp", record.id(), record.seq_len());
//!
//! // Fetch a specific region
//! let region = reader.fetch_range("chr1", 1000, 2000).unwrap();
//! println!("Region: {} bp", region.len());
//! ```

// Indexed random access modules
pub mod bgzf;
pub mod fai;
pub mod gzi;
pub mod indexed;

#[cfg(feature = "url")]
pub mod remote;

#[allow(non_snake_case)]
pub mod FastX
{
    use flate2::read::MultiGzDecoder;
    use std::ffi::OsStr;
    use std::io;
    use std::io::BufRead;

    const PER_THREAD_BUF_SIZE: usize = 600 * 1024 * 1024;

    /// Sequence file format type.
    ///
    /// Represents the detected format of a sequence file based on its first byte.
    pub enum FastXFormat
    {
        /// FASTQ format - records starting with `@`
        FASTQ,
        /// FASTA format - records starting with `>`
        FASTA,
        /// End of file reached
        EOF,
        /// Unknown or invalid format
        UNKNOWN,
    }

    /// A FASTA sequence record.
    ///
    /// FASTA format is a text-based format for representing nucleotide or peptide sequences.
    /// Each record consists of a header line starting with `>` followed by sequence data
    /// that may span multiple lines.
    ///
    /// # Fields
    ///
    /// * `name` - The full header line (without the leading `>`)
    /// * `raw_seq` - The raw sequence data including newlines (for multi-line sequences)
    ///
    /// # Example
    ///
    /// ```
    /// use fastx::FastX::{FastARecord, FastXRead};
    /// use std::io::BufReader;
    ///
    /// let data = b">seq1 description\nACGT\nACGT\n";
    /// let mut reader = BufReader::new(&data[..]);
    /// let mut record = FastARecord::default();
    /// record.read(&mut reader).unwrap();
    ///
    /// assert_eq!(record.name(), "seq1 description");
    /// assert_eq!(record.id(), "seq1");
    /// assert_eq!(record.seq(), b"ACGTACGT");
    /// ```
    #[derive(Default)]
    pub struct FastARecord
    {
        /// Full header line (without leading `>`)
        pub name: String,
        /// Raw sequence data (may contain newlines for multi-line sequences)
        pub raw_seq: Vec<u8>,
    }

    /// A FASTQ sequence record.
    ///
    /// FASTQ format extends FASTA by adding base quality scores. Each record has four lines:
    /// 1. Header line starting with `@`
    /// 2. Sequence data
    /// 3. Separator line starting with `+` (optionally repeating the header)
    /// 4. Quality scores (same length as sequence)
    ///
    /// Quality scores are encoded as ASCII characters, with each character representing
    /// the quality of the corresponding base in the sequence.
    ///
    /// # Example
    ///
    /// ```
    /// use fastx::FastX::{FastQRecord, FastXRead, FastQRead};
    /// use std::io::BufReader;
    ///
    /// let data = b"@seq1\nACGT\n+\n!!!!\n";
    /// let mut reader = BufReader::new(&data[..]);
    /// let mut record = FastQRecord::default();
    /// record.read(&mut reader).unwrap();
    ///
    /// assert_eq!(record.name(), "seq1");
    /// assert_eq!(record.seq(), b"ACGT");
    /// assert_eq!(record.qual(), &b"!!!!".to_vec());
    /// ```
    #[derive(Default)]
    pub struct FastQRecord
    {
        name: String,
        seq: Vec<u8>,
        comment: String,
        qual: Vec<u8>,
    }

    /// Core trait for reading FASTA/FASTQ records.
    ///
    /// This trait provides methods for reading sequence records and accessing their components.
    /// It is implemented by both [`FastARecord`] and [`FastQRecord`].
    pub trait FastXRead: std::fmt::Display
    {
        /// Read the next record from the reader.
        ///
        /// # Returns
        ///
        /// * `Ok(n)` where `n` is the number of bytes read. Returns `0` on EOF.
        /// * `Err(e)` if an I/O error occurred or the data is malformed.
        fn read(&mut self, reader: &mut dyn BufRead) -> io::Result<usize>;

        /// Get the full header/name line (without the leading `>` or `@`).
        fn name(&self) -> &String;

        /// Get the sequence identifier (the part before the first space in the name).
        ///
        /// For FASTA, this is the name without the leading `>`.
        /// For FASTQ, this is the name without the leading `@`.
        fn id(&self) -> &str;

        /// Get the description (the part after the first space in the name).
        ///
        /// Returns an empty string if there is no description.
        fn desc(&self) -> &str;

        /// Get the raw sequence data including any newlines.
        ///
        /// For FASTA records with multi-line sequences, this preserves the newlines.
        /// For FASTQ records, this is typically a single line.
        fn seq_raw(&self) -> &Vec<u8>;

        /// Get the sequence data with newlines removed.
        ///
        /// This returns a contiguous sequence by stripping out line breaks.
        fn seq(&self) -> Vec<u8>;

        /// Get the sequence length (excluding newlines).
        fn seq_len(&self) -> usize;

        /// Get the sequence split into individual lines.
        fn lines(&self) -> Vec<&[u8]>;
    }

    /// Iterator wrapper for reading sequence records.
    ///
    /// This provides an iterator interface over FASTA or FASTQ records.
    /// Use [`fasta_iter`] or [`fastq_iter`] to create instances.
    ///
    /// # Type Parameters
    ///
    /// * `R` - The buffered reader type
    /// * `T` - The record type ([`FastARecord`] or [`FastQRecord`])
    ///
    /// # Example
    ///
    /// ```no_run
    /// use fastx::FastX::{fasta_iter, FastXRead};
    /// use std::io::BufReader;
    /// use std::fs::File;
    ///
    /// let reader = BufReader::new(File::open("sequences.fasta").unwrap());
    /// for result in fasta_iter(reader) {
    ///     let record = result.unwrap();
    ///     println!("Sequence: {}", record.id());
    /// }
    /// ```
    pub struct FastXIterator<R: BufRead, T: FastXRead + Default>
    {
        reader: R,
        done: bool,
        _phantom: std::marker::PhantomData<T>,
    }

    impl<R: BufRead, T: FastXRead + Default> FastXIterator<R, T>
    {
        /// Create a new iterator over sequence records.
        ///
        /// Generally, you should use [`fasta_iter`] or [`fastq_iter`] instead
        /// of calling this directly.
        pub fn new(reader: R) -> Self
        {
            Self {
                reader,
                done: false,
                _phantom: std::marker::PhantomData,
            }
        }
    }

    impl<R: BufRead, T: FastXRead + Default> Iterator for FastXIterator<R, T>
    {
        type Item = io::Result<T>;

        fn next(&mut self) -> Option<Self::Item>
        {
            if self.done
            {
                return None;
            }

            let mut new_record = T::default();
            match new_record.read(&mut self.reader)
            {
                Ok(0) =>
                {
                    self.done = true;
                    None
                }
                Ok(_) => Some(Ok(new_record)),
                Err(e) =>
                {
                    self.done = true;
                    Some(Err(e))
                }
            }
        }
    }

    /// Trait for reading FASTQ-specific data.
    ///
    /// This trait extends [`FastXRead`] with methods specific to FASTQ format,
    /// which includes quality scores for each base in the sequence.
    ///
    /// Only [`FastQRecord`] implements this trait.
    pub trait FastQRead: FastXRead
    {
        /// Get the comment line (the content after the `+` separator).
        ///
        /// This is typically empty, but some FASTQ files repeat the sequence name here.
        fn comment(&self) -> &str;

        /// Get the quality scores.
        ///
        /// Returns a byte vector where each byte represents the quality score
        /// for the corresponding base in the sequence. Quality scores are typically
        /// encoded in Phred format.
        fn qual(&self) -> &Vec<u8>;
    }

    impl std::fmt::Display for FastARecord
    {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result
        {
            write!(f, ">{}\n{}", self.name(), String::from_utf8_lossy(&self.seq()))
        }
    }

    impl FastXRead for FastARecord
    {
        fn name(&self) -> &String
        {
            &self.name
        }

        fn id(&self) -> &str
        {
            match memchr::memchr(b' ', self.name.as_bytes())
            {
                None => &self.name,
                Some(i) => &self.name[..i],
            }
        }

        fn desc(&self) -> &str
        {
            match memchr::memchr(b' ', self.name.as_bytes())
            {
                None => "",
                Some(i) => &self.name[i + 1..],
            }
        }

        fn seq_raw(&self) -> &Vec<u8>
        {
            &self.raw_seq
        }

        fn seq(&self) -> Vec<u8>
        {
            let mut seq = vec![0; self.raw_seq.len()];
            let mut line_start = 0;
            let mut seq_end = 0;
            let mut seq_start = 0;
            memchr::memchr_iter(b'\n', &self.raw_seq).for_each(|line_end| {
                seq_start = seq_end;
                seq_end += line_end - line_start;
                seq[seq_start..seq_end].copy_from_slice(&self.raw_seq[line_start..line_end]);
                line_start = line_end + 1; // skip '\n'
            });
            if line_start < self.raw_seq.len()
            {
                seq_start = seq_end;
                seq_end += self.raw_seq.len() - line_start;
                seq[seq_start..seq_end]
                    .copy_from_slice(&self.raw_seq[line_start..self.raw_seq.len()]);
                seq.resize(seq_end, 0);
            }
            seq
        }

        fn lines(&self) -> Vec<&[u8]>
        {
            //self.seq.split( |c| *c == b'\n').collect()
            let mut line_start = 0;
            memchr::memchr_iter(b'\n', &self.raw_seq)
                .map(|line_end| {
                    let line = &self.raw_seq[line_start..line_end];
                    line_start = line_end + 1;
                    line
                })
                .collect()
        }

        fn seq_len(&self) -> usize
        {
            let mut line_start = 0;
            memchr::memchr_iter(b'\n', &self.raw_seq).fold(0, |mut len, line_end| {
                len += line_end - line_start;
                line_start = line_end + 1;
                len
            }) + self.raw_seq.len()
                - line_start
        }

        fn read(&mut self, reader: &mut dyn BufRead) -> io::Result<usize>
        {
            self.name.clear();
            self.raw_seq.clear();

            let mut size = 0;

            let mut record_sep = [0_u8];
            match reader.read(&mut record_sep)
            {
                Err(e) => return Err(e),
                Ok(0) => return Ok(0),
                Ok(some) =>
                {
                    size += some;
                    if some != 1 && record_sep[0] != b'>'
                    {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "record_sep does not match '>'",
                        ));
                    }
                }
            };

            match reader.read_line(&mut self.name)
            {
                Err(e) => return Err(e),
                Ok(0) => return Ok(0),
                Ok(some) => size += some,
            };
            rstrip_newline_string(&mut self.name);

            match read_until_before(reader, b'>', &mut self.raw_seq)
            {
                Err(e) => Err(e),
                Ok(0) => Ok(0),
                Ok(some) =>
                {
                    rstrip_seq(&mut self.raw_seq);
                    Ok(size + some)
                }
            }
        }
    }

    impl std::fmt::Display for FastQRecord
    {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result
        {
            write!(
                f,
                "@{}\n{}\n+\n{}",
                self.name(),
                String::from_utf8_lossy(&self.seq()),
                String::from_utf8_lossy(&self.qual())
            )
        }
    }

    impl FastXRead for FastQRecord
    {
        fn name(&self) -> &String
        {
            &self.name
        }

        fn id(&self) -> &str
        {
            match memchr::memchr(b' ', self.name.as_bytes())
            {
                None => &self.name,
                Some(i) => &self.name[..i],
            }
        }

        fn desc(&self) -> &str
        {
            match memchr::memchr(b' ', self.name.as_bytes())
            {
                None => "",
                Some(i) => &self.name[i + 1..],
            }
        }

        fn seq_raw(&self) -> &Vec<u8>
        {
            &self.seq
        }

        // As multiline FastQ is very uncommon, we assume seq to be one line
        fn seq(&self) -> Vec<u8>
        {
            self.seq.clone()
        }

        fn seq_len(&self) -> usize
        {
            self.seq
                .split(|c| *c == b'\n')
                .fold(0, |len, seq| len + seq.len())
        }

        fn lines(&self) -> Vec<&[u8]>
        {
            self.seq.split(|c| *c == b'\n').collect()
        }

        fn read(&mut self, reader: &mut dyn BufRead) -> io::Result<usize>
        {
            self.name.clear();
            let mut size;
            match reader.read_line(&mut self.name)
            {
                Err(e) => return Err(e),
                Ok(0) => return Ok(0),
                Ok(some) => size = some,
            }
            rstrip_newline_string(&mut self.name);

            if self.name.is_empty() {
                 return Err(io::Error::new(io::ErrorKind::InvalidData, "FASTQ record with empty header"));
            }

            if self.name.starts_with('@') {
                self.name.remove(0);
            } else {
                 return Err(io::Error::new(io::ErrorKind::InvalidData, "FASTQ header must start with @"));
            }

            self.seq.clear();
            match reader.read_until(b'\n', &mut self.seq)
            {
                Err(e) => return Err(e),
                Ok(0) => return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "FASTQ truncated sequence")),
                Ok(some) =>
                {
                    rstrip_newline_vec(&mut self.seq);
                    size += some;
                }
            }

            self.comment.clear();
            match reader.read_line(&mut self.comment)
            {
                Err(e) => return Err(e),
                Ok(0) => return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "FASTQ truncated comment")),
                Ok(some) =>
                {
                    rstrip_newline_string(&mut self.comment);
                    if !self.comment.starts_with('+') {
                        return Err(io::Error::new(io::ErrorKind::InvalidData, "FASTQ comment must start with +"));
                    }
                    size += some
                }
            }

            self.qual.clear();
            match reader.read_until(b'\n', &mut self.qual)
            {
                Err(e) => Err(e),
                Ok(0) => Err(io::Error::new(io::ErrorKind::UnexpectedEof, "FASTQ truncated quality")),
                Ok(some) =>
                {
                    rstrip_newline_vec(&mut self.qual);
                    Ok(size + some)
                }
            }
        }
    }

    impl FastQRead for FastQRecord
    {
        fn comment(&self) -> &str
        {
            if self.comment.len() > 0 {
                &self.comment[1..]
            } else {
                ""
            }
        }

        fn qual(&self) -> &Vec<u8>
        {
            &self.qual
        }
    }

    fn rstrip_newline_string(s: &mut String)
    {
        while s.ends_with('\n') || s.ends_with('\r')
        {
            s.pop();
        }
    }

    fn rstrip_seq(s: &mut Vec<u8>)
    {
        while s.last() == Some(&b'>') || s.last() == Some(&b'\n') || s.last() == Some(&b'\r')
        {
            s.pop();
        }
    }

    fn rstrip_newline_vec(s: &mut Vec<u8>)
    {
        while s.last() == Some(&b'\n') || s.last() == Some(&b'\r')
        {
            s.pop();
        }
    }

    /// Peek at the first byte to determine the file format.
    ///
    /// This function reads (but does not consume) the first byte of the input
    /// to determine whether it's a FASTA or FASTQ file.
    ///
    /// # Arguments
    ///
    /// * `reader` - A buffered reader positioned at the start of a sequence file
    ///
    /// # Returns
    ///
    /// * `Ok((format, first_byte))` - The detected format and the first byte
    /// * `Err(e)` - If the format is unknown or invalid
    ///
    /// # Format Detection
    ///
    /// * `>` - FASTA format
    /// * `@` - FASTQ format
    /// * `\0` - EOF
    /// * Other - UNKNOWN (returns error)
    ///
    /// # Example
    ///
    /// ```
    /// use fastx::FastX::{peek, FastXFormat};
    /// use std::io::BufReader;
    ///
    /// let data = b">sequence\nACGT\n";
    /// let mut reader = BufReader::new(&data[..]);
    /// let (format, _) = peek(&mut reader).unwrap();
    /// matches!(format, FastXFormat::FASTA);
    /// ```
    pub fn peek(reader: &mut dyn BufRead) -> io::Result<(FastXFormat, u8)>
    {
        let buf = reader.fill_buf().expect("peek failed");
        let format = match buf[0] as char
        {
            '>' => FastXFormat::FASTA,
            '@' => FastXFormat::FASTQ,
            '\0' => FastXFormat::EOF,
            _ => FastXFormat::UNKNOWN,
        };
        if let FastXFormat::UNKNOWN = format
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Wrong format expected '>' or '@'!",
            ));
        }
        Ok((format, buf[0]))
    }

    use std::fs::File;
    use std::io::BufReader;
    use std::path::Path;
    //use std::str::pattern::Pattern;

    /// Create a buffered reader from a file path.
    ///
    /// This function opens a file and wraps it in a large (600MB) buffered reader
    /// for optimal performance. If the file has a `.gz` extension, it automatically
    /// applies gzip decompression.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the sequence file (can be `.fasta`, `.fastq`, `.fasta.gz`, `.fastq.gz`, etc.)
    ///
    /// # Returns
    ///
    /// A boxed buffered reader ready for reading sequence data.
    ///
    /// # Example
    ///
    /// ```no_run
    /// use fastx::FastX::{reader_from_path, fasta_iter, FastXRead};
    /// use std::path::Path;
    ///
    /// // Automatically handles both compressed and uncompressed files
    /// let reader = reader_from_path(Path::new("sequences.fasta.gz")).unwrap();
    /// for record in fasta_iter(reader) {
    ///     println!("{}", record.unwrap().id());
    /// }
    /// ```
    pub fn reader_from_path(path: &Path) -> io::Result<Box<dyn BufRead>>
    {
        let file = File::open(path)?;
        let reader: Box<dyn BufRead> = match path.extension()
        {
            Some(extension) if extension == OsStr::new("gz") => Box::new(BufReader::with_capacity(
                PER_THREAD_BUF_SIZE,
                MultiGzDecoder::new(BufReader::new(file)),
            )),
            _ => Box::new(BufReader::with_capacity(PER_THREAD_BUF_SIZE, file)),
        };
        Ok(reader)
    }

    /// Create a record reader with automatic format detection.
    ///
    /// This function peeks at the first byte of the input to determine whether
    /// it's FASTA or FASTQ format, then returns an appropriate record reader.
    ///
    /// # Arguments
    ///
    /// * `reader` - A buffered reader positioned at the start of sequence data
    ///
    /// # Returns
    ///
    /// A boxed trait object implementing [`FastXRead`] (either [`FastARecord`] or [`FastQRecord`]).
    ///
    /// # Example
    ///
    /// ```no_run
    /// use fastx::FastX::{from_reader, FastXRead};
    /// use std::io::BufReader;
    /// use std::fs::File;
    ///
    /// let file = File::open("unknown_format.txt").unwrap();
    /// let mut reader = BufReader::new(file);
    /// let mut record = from_reader(&mut reader).unwrap();
    ///
    /// while record.read(&mut reader).unwrap() > 0 {
    ///     println!("{}", record.id());
    /// }
    /// ```
    pub fn from_reader(reader: &mut dyn BufRead) -> io::Result<Box<dyn FastXRead>>
    {
        let (format, first) = peek(reader)?;
        match format
        {
            FastXFormat::FASTA => Ok(Box::new(FastARecord::default())),
            FastXFormat::FASTQ => Ok(Box::new(FastQRecord::default())),
            FastXFormat::EOF | FastXFormat::UNKNOWN =>
            {
                Err(io::Error::new(io::ErrorKind::InvalidData, format!("{:?}", first)))
            }
        }
    }

    /// Create an iterator over FASTA records.
    ///
    /// This is the recommended way to read FASTA files. Each iteration yields
    /// a `Result<FastARecord, io::Error>`.
    ///
    /// # Arguments
    ///
    /// * `reader` - A buffered reader containing FASTA data
    ///
    /// # Returns
    ///
    /// A [`FastXIterator`] that yields FASTA records.
    ///
    /// # Example
    ///
    /// ```no_run
    /// use fastx::FastX::{fasta_iter, FastXRead};
    /// use std::io::BufReader;
    /// use std::fs::File;
    ///
    /// let reader = BufReader::new(File::open("sequences.fasta").unwrap());
    /// for result in fasta_iter(reader) {
    ///     let record = result.unwrap();
    ///     println!("{}: {} bp", record.id(), record.seq_len());
    /// }
    /// ```
    pub fn fasta_iter<R: BufRead>(reader: R) -> FastXIterator<R, FastARecord>
    {
        FastXIterator::new(reader)
    }

    /// Create an iterator over FASTQ records.
    ///
    /// This is the recommended way to read FASTQ files. Each iteration yields
    /// a `Result<FastQRecord, io::Error>`.
    ///
    /// # Arguments
    ///
    /// * `reader` - A buffered reader containing FASTQ data
    ///
    /// # Returns
    ///
    /// A [`FastXIterator`] that yields FASTQ records.
    ///
    /// # Example
    ///
    /// ```no_run
    /// use fastx::FastX::{fastq_iter, FastXRead, FastQRead};
    /// use std::io::BufReader;
    /// use std::fs::File;
    ///
    /// let reader = BufReader::new(File::open("sequences.fastq").unwrap());
    /// for result in fastq_iter(reader) {
    ///     let record = result.unwrap();
    ///     println!("{}: {} bp, {} qual", record.id(), record.seq_len(), record.qual().len());
    /// }
    /// ```
    pub fn fastq_iter<R: BufRead>(reader: R) -> FastXIterator<R, FastQRecord>
    {
        FastXIterator::new(reader)
    }

    /// Iterate over FASTA records with buffer reuse for high performance.
    ///
    /// This function calls the provided closure for each record, reusing the same
    /// buffer memory to avoid allocations. This is significantly faster than
    /// [`fasta_iter`] but does not allow the records to outlive the closure.
    ///
    /// # Arguments
    ///
    /// * `reader` - A buffered reader containing FASTA data
    /// * `func` - A closure that takes a reference to a [`FastARecord`]
    ///
    /// # Example
    ///
    /// ```no_run
    /// use fastx::FastX::{fasta_for_each, FastXRead};
    /// use std::io::BufReader;
    /// use std::fs::File;
    ///
    /// let reader = BufReader::new(File::open("sequences.fasta").unwrap());
    /// fasta_for_each(reader, |record| {
    ///     println!("{}: {} bp", record.id(), record.seq_len());
    /// }).unwrap();
    /// ```
    pub fn fasta_for_each<R: BufRead, F>(mut reader: R, mut func: F) -> io::Result<()>
    where
        F: FnMut(&FastARecord),
    {
        let mut record = FastARecord::default();
        while record.read(&mut reader)? > 0
        {
            func(&record);
        }
        Ok(())
    }

    /// Iterate over FASTQ records with buffer reuse for high performance.
    ///
    /// This function calls the provided closure for each record, reusing the same
    /// buffer memory to avoid allocations. This is significantly faster than
    /// [`fastq_iter`] but does not allow the records to outlive the closure.
    ///
    /// # Arguments
    ///
    /// * `reader` - A buffered reader containing FASTQ data
    /// * `func` - A closure that takes a reference to a [`FastQRecord`]
    ///
    /// # Example
    ///
    /// ```no_run
    /// use fastx::FastX::{fastq_for_each, FastXRead};
    /// use std::io::BufReader;
    /// use std::fs::File;
    ///
    /// let reader = BufReader::new(File::open("sequences.fastq").unwrap());
    /// fastq_for_each(reader, |record| {
    ///     println!("{}: {} bp", record.id(), record.seq_len());
    /// }).unwrap();
    /// ```
    pub fn fastq_for_each<R: BufRead, F>(mut reader: R, mut func: F) -> io::Result<()>
    where
        F: FnMut(&FastQRecord),
    {
        let mut record = FastQRecord::default();
        while record.read(&mut reader)? > 0
        {
            func(&record);
        }
        Ok(())
    }

    /// Iterate over sequence records with automatic format detection and buffer reuse.
    ///
    /// This function peeks at the first byte to determine if the file is FASTA or FASTQ,
    /// then calls the appropriate `for_each` function.
    ///
    /// # Arguments
    ///
    /// * `reader` - A buffered reader containing FASTA or FASTQ data
    /// * `fasta_func` - A closure that takes a reference to a [`FastARecord`]
    /// * `fastq_func` - A closure that takes a reference to a [`FastQRecord`]
    ///
    /// # Example
    ///
    /// ```no_run
    /// use fastx::FastX::{fastx_for_each, FastXRead};
    /// use std::io::BufReader;
    /// use std::fs::File;
    ///
    /// let mut reader = BufReader::new(File::open("sequences.fasta").unwrap());
    /// fastx_for_each(reader,
    ///     |record| println!("FASTA: {}", record.id()),
    ///     |record| println!("FASTQ: {}", record.id()),
    /// ).unwrap();
    /// ```
    pub fn fastx_for_each<R: BufRead, FA, FQ>(
        mut reader: R,
        fasta_func: FA,
        fastq_func: FQ,
    ) -> io::Result<()>
    where
        FA: FnMut(&FastARecord),
        FQ: FnMut(&FastQRecord),
    {
        let (format, _) = peek(&mut reader)?;
        match format
        {
            FastXFormat::FASTA => fasta_for_each(reader, fasta_func),
            FastXFormat::FASTQ => fastq_for_each(reader, fastq_func),
            FastXFormat::EOF | FastXFormat::UNKNOWN => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Unknown sequence format",
            )),
        }
    }

    /// from std::io::read_until, adapted to not consume the delimiter
    fn read_until_before<R: BufRead + ?Sized>(
        r: &mut R,
        delim: u8,
        buf: &mut Vec<u8>,
    ) -> io::Result<usize>
    {
        let mut read = 0;
        loop
        {
            let (done, used) = {
                let available = match r.fill_buf()
                {
                    Ok(n) => n,
                    //Err(ref e) if e.is_interrupted() => continue,
                    Err(e) => return Err(e),
                };
                match memchr::memchr(delim, available)
                {
                    Some(i) =>
                    {
                        buf.extend_from_slice(&available[..=i]);
                        (true, i + 1)
                    }
                    None =>
                    {
                        buf.extend_from_slice(available);
                        (false, available.len())
                    }
                }
            };
            if done
            {
                r.consume(used - 1); // do not consume delimiter
                read += used - 1;
                return Ok(read);
            }
            else
            {
                r.consume(used);
                read += used;
            }
            if used == 0
            {
                return Ok(read);
            }
        }
    }
}

#[cfg(test)]
mod tests
{
    use super::FastX::FastARecord;
    use super::FastX::FastQRecord;
    use super::FastX::FastXRead;
    use std::io::BufReader;
    use std::io::Cursor;

    #[test]
    fn fasta()
    {
        let mut x = BufReader::new(Cursor::new(">a\nAGTC\n>b\nTAGC\nTTTT\n>c\nGCTA"));
        let mut record = FastARecord::default();
        let _ = record.read(&mut x);
        assert_eq!("a", record.name());
        assert_eq!(b"AGTC".to_vec(), record.seq());
        assert_eq!(&b"AGTC".to_vec(), record.seq_raw());
        let _ = record.read(&mut x);
        assert_eq!("b", record.name());
        assert_eq!(b"TAGCTTTT".to_vec(), record.seq());
        assert_eq!(&b"TAGC\nTTTT".to_vec(), record.seq_raw());
        let _ = record.read(&mut x);
        assert_eq!("c", record.name());
        assert_eq!(b"GCTA".to_vec(), record.seq());
        assert_eq!(&b"GCTA".to_vec(), record.seq_raw());
    }

    #[test]
    fn fasta_memchr()
    {
        let fasta = b">a\nAGTC\n>b\nTAGC\nTTTT\n>c\nGCTA\n";
        let mut reader = BufReader::new(Cursor::new(fasta));
        let mut record = FastARecord::default();
        let mut output: Vec<u8> = Vec::new();
        //peek(&mut reader).expect("peek");
        while let Ok(_some @ 1..=usize::MAX) = record.read(&mut reader)
        {
            output.push(b'>');
            output.append(&mut record.name().as_bytes().to_vec());
            output.push(b'\n');
            let seq = record.seq_raw();
            println!("r{}r", String::from_utf8_lossy(seq));
            let mut offset = 0;
            memchr::memchr_iter(b'\n', seq)
                .map(|line_end| {
                    let seq = &seq[offset..line_end];
                    println!("#{}#", String::from_utf8_lossy(seq));
                    offset = line_end + 1;
                    seq
                })
                .for_each(|line| {
                    output.append(&mut line.to_vec());
                    output.push(b'\n');
                });
            output.append(&mut seq[offset..].to_vec());
            output.push(b'\n');
        }
        println!("{}\n\n{}", String::from_utf8_lossy(fasta), String::from_utf8_lossy(&output));
        assert_eq!(fasta.to_vec(), output);
    }

    #[test]
    fn fastq()
    {
        let mut x = BufReader::new(Cursor::new(
            "@a\nAGTC\n+\n'&'*+\n@b\nTAGCTTTT\n+\n'&'*+'&'*+\n@c\nGCTA\n+\n'&'*+",
        ));
        let mut record = FastQRecord::default();
        let _ = record.read(&mut x);
        assert_eq!("a", record.name());
        assert_eq!(b"AGTC".to_vec(), record.seq());
        assert_eq!(&b"AGTC".to_vec(), record.seq_raw());

        let _ = record.read(&mut x);
        assert_eq!("b", record.name());
        assert_eq!(b"TAGCTTTT".to_vec(), record.seq());
        assert_eq!(&b"TAGCTTTT".to_vec(), record.seq_raw());

        let _ = record.read(&mut x);
        assert_eq!("c", record.name());
        assert_eq!(b"GCTA".to_vec(), record.seq());
        assert_eq!(&b"GCTA".to_vec(), record.seq_raw());
    }

    #[test]
    fn fasta_iterator()
    {
        use super::FastX::fasta_iter;
        let reader = BufReader::new(Cursor::new(">a\nAGTC\n>b\nTAGC\nTTTT\n>c\nGCTA"));
        let records: Result<Vec<_>, _> = fasta_iter(reader).collect();
        let records = records.unwrap();

        assert_eq!(3, records.len());
        assert_eq!("a", records[0].name());
        assert_eq!(b"AGTC".to_vec(), records[0].seq());
        assert_eq!("b", records[1].name());
        assert_eq!(b"TAGCTTTT".to_vec(), records[1].seq());
        assert_eq!("c", records[2].name());
        assert_eq!(b"GCTA".to_vec(), records[2].seq());
    }

    #[test]
    fn fastq_iterator()
    {
        use super::FastX::fastq_iter;
        let reader = BufReader::new(Cursor::new(
            "@a\nAGTC\n+\n'&'*+\n@b\nTAGCTTTT\n+\n'&'*+'&'*+\n@c\nGCTA\n+\n'&'*+",
        ));
        let records: Result<Vec<_>, _> = fastq_iter(reader).collect();
        let records = records.unwrap();

        assert_eq!(3, records.len());
        assert_eq!("a", records[0].name());
        assert_eq!(b"AGTC".to_vec(), records[0].seq());
        assert_eq!("b", records[1].name());
        assert_eq!(b"TAGCTTTT".to_vec(), records[1].seq());
        assert_eq!("c", records[2].name());
        assert_eq!(b"GCTA".to_vec(), records[2].seq());
    }
}
