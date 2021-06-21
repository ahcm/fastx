#[allow(non_snake_case)]
pub mod FastX {
    use std::io;
    use std::io::BufRead;
    use std::io::Read;

    pub enum FastXFormat
    {
        FASTQ,
        FASTA,
        UNKNOWN,
    }

    #[derive(Default)]
    pub struct FastXRecord
    {
        pub name: String,
        pub seq: Vec<u8>,
    }

    #[derive(Default)]
    pub struct FastQRecord
    {
        name: String,
        seq: Vec<u8>,
        comment: String,
        qual: Vec<u8>,
    }

    pub trait FastXRead
    {
        fn read(&mut self, reader: &mut dyn BufRead) -> io::Result<usize>;
        fn name(&self) -> &String;
        fn seq_raw(&self) -> &Vec<u8>;
        fn seq(&self) -> Vec<u8>;
        fn seq_len(&self) -> usize;
        fn lines(&self) -> Vec<&[u8]>;
    }

    pub trait FastQRead: FastXRead
    {
        fn qual(&self) -> &Vec<u8>;
    }

    impl FastXRead for FastXRecord
    {
        fn name(&self) -> &String
        {
            &self.name
        }

        fn seq_raw(&self) -> &Vec<u8>
        {
            &self.seq
        }

        fn seq(&self) -> Vec<u8>
        {
            /*
            self.seq
                .lines()
                .fold(String::with_capacity(self.seq.len()), |mut acc, line| {
                    acc.push_str(&line.unwrap());
                    acc
                })
                */
            
            let mut seq = Vec::with_capacity(self.seq.len());
            seq.resize(self.seq.len(), 0);
            let mut line_start = 0;
            let mut seq_end = 0;
            let mut seq_start = 0;
            memchr::memchr_iter(b'\n', &self.seq)
                .for_each(|line_end|
                     {
                         seq_start = seq_end;
                         seq_end += line_end - line_start;
                         seq[seq_start..seq_end].copy_from_slice(&self.seq[line_start..line_end]);
                         line_start = line_end + 1; // skip '\n'
                     });
            if line_start < self.seq.len()
            {
                seq_start = seq_end;
                seq_end += self.seq.len() - line_start;
                seq[seq_start..seq_end].copy_from_slice(&self.seq[line_start..self.seq.len()]);
                seq.resize(seq_end, 0);
            }
            seq
        }

        fn lines(&self) -> Vec<&[u8]>
        {
            //self.seq.split( |c| *c == b'\n').collect()
            let mut line_start = 0;
            memchr::memchr_iter(b'\n', &self.seq)
                .map(|line_end|
                     {
                         let line = &self.seq[line_start..line_end];
                         line_start = line_end + 1;
                         line
                     }).collect()
        }

        fn seq_len(&self) -> usize
        {
            let mut line_start = 0;
            memchr::memchr_iter(b'\n', &self.seq)
                .fold(0, |mut len, line_end|
                     {
                         len += line_end - line_start;
                         line_start = line_end + 1;
                         len
                     })
                + self.seq.len() - line_start
        }

        fn read(&mut self, reader: &mut dyn BufRead) -> io::Result<usize> {
            self.name.clear();
            let size;
            match reader.read_line(&mut self.name) {
                Err(e) => return Err(e),
                Ok(0) => return Ok(0),
                Ok(some) => size = some,
            }
            rstrip_newline_string(&mut self.name);

            self.seq.clear();
            match reader.read_until(b'>', &mut self.seq) {
                Err(e) => Err(e),
                Ok(0) => Ok(0),
                Ok(some) => {
                    rstrip_seq(&mut self.seq);
                    //println!("{}", String::from_utf8_lossy(&self.seq));
                    Ok(size + some)
                }
            }
        }
    }

    impl FastXRead for FastQRecord {
        fn name(&self) -> &String
        {
            &self.name
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
            self.seq.split( |c| *c == b'\n').fold(0, |len, seq| len + seq.len())
        }

        fn lines(&self) -> Vec<&[u8]>
        {
            self.seq.split( |c| *c == b'\n').collect()
        }

        fn read(&mut self, reader: &mut dyn BufRead) -> io::Result<usize> {
            self.name.clear();
            let mut size;
            match reader.read_line(&mut self.name) {
                Err(e) => return Err(e),
                Ok(0) => return Ok(0),
                Ok(some) => size = some,
            }
            rstrip_newline_string(&mut self.name); //self.name.truncate(size - 1); // truncate newline XXX non UNIX

            self.seq.clear();
            match reader.read_until(b'\n', &mut self.seq) {
                Err(e) => return Err(e),
                Ok(0) => return Ok(0),
                Ok(some) => {
                    rstrip_newline_vec(&mut self.seq);
                    //println!("{}", String::from_utf8_lossy(&self.seq));
                    size += some;
                }
            }

            self.comment.clear();
            match reader.read_line(&mut self.comment) {
                Err(e) => return Err(e),
                Ok(0) => return Ok(0),
                Ok(some) => size += some,
            }
            rstrip_newline_string(&mut self.comment); //self.name.truncate(size - 1); // truncate newline XXX non UNIX

            self.qual.clear();
            match reader.read_until(b'\n', &mut self.qual) {
                Err(e) => Err(e),
                Ok(0) => Ok(0),
                Ok(some) => {
                    rstrip_newline_vec(&mut self.seq);
                    //println!("{}", String::from_utf8_lossy(&self.seq));
                    Ok(size + some)
                }
            }
        }
    }

    impl FastQRead for FastQRecord {
        fn qual(&self) -> &Vec<u8> {
            &self.qual
        }
    }

    fn rstrip_newline_string(s: &mut String) {
        while s.ends_with(&"\n") {
            s.truncate(s.len() - 1);
        }
    }

    fn rstrip_seq(s: &mut Vec<u8>) {
        while s[s.len() - 1] == b'>' || s[s.len() - 1] == b'\n'
        //.ends_with(&[b'>',b'\n'])
        {
            s.truncate(s.len() - 1);
        }
    }

    fn rstrip_newline_vec(s: &mut Vec<u8>) {
        while s[s.len() - 1] == b'\n' {
            s.truncate(s.len() - 1);
        }
    }

    // Determines the format from reading the first byte(s) of
    // the Readable
    pub fn peek(reader: &mut dyn Read) -> io::Result<(FastXFormat, [u8; 1])> {
        let mut buf = [0 as u8];
        reader.read_exact(&mut buf)?;
        let format = match buf[0] as char {
            '>' => FastXFormat::FASTA,
            '@' => FastXFormat::FASTQ,
            _ => FastXFormat::UNKNOWN,
        };
        if let FastXFormat::UNKNOWN = format {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Wrong format expected '>' or '@'!",
            ));
        }
        Ok((format, buf))
    }

    pub fn open(reader: &mut dyn Read) -> io::Result<Box<dyn FastXRead>> {
        let (format, first) = peek(reader)?;
        match format {
            FastXFormat::FASTA => Ok(Box::new(FastXRecord::default())),
            FastXFormat::FASTQ => Ok(Box::new(FastQRecord::default())),
            FastXFormat::UNKNOWN => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("{:?}", first),
            )),
        }
    }

    //            let mut seq = String::new();
    //            let mut next_line = String::new();
    //            loop
    //            {
    //
    //                match reader.read_until() read_line(&mut next_line)
    //                {
    //                    Err(e) => return Err(e),
    //                    Ok(0) => break,
    //                    Ok(_some) =>
    //                        if next_line.starts_with('>')
    //                        {
    //                            break;
    //                        }
    //                    else
    //                    {
    //                        seq.extend()
    //                    }
    //                }
    //
    //            }
    //        }
    //    }
}

#[cfg(test)]
mod tests {
    use super::FastX::FastQRecord;
    use super::FastX::FastXRead;
    use super::FastX::FastXRecord;
    use super::FastX::peek;
    use std::io::BufReader;
    use std::io::Cursor;

    #[test]
    fn fasta() {
        let mut x = BufReader::new(Cursor::new(">a\nAGTC\n>b\nTAGC\nTTTT\n>c\nGCTA"));
        let mut record = FastXRecord::default();
        let _ = record.read(&mut x);
        assert_eq!(">a", record.name());
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
        let mut record = FastXRecord::default();
        let mut output : Vec<u8> = Vec::new();
        peek(&mut reader).expect("peek");
        while let Ok(_some @ 1..=usize::MAX) = record.read(&mut reader)
        {
            output.push(b'>');
            output.append(&mut record.name().as_bytes().to_vec());
            output.push(b'\n');
            let seq = record.seq_raw();
            println!("r{}r", String::from_utf8_lossy(seq));
            let mut offset = 0;
            memchr::memchr_iter(b'\n', &seq)
                .map(|line_end|
                     {
                         let seq = &seq[offset..line_end];
                         println!("#{}#", String::from_utf8_lossy(seq));
                         offset = line_end + 1;
                         seq
                     })
                .for_each(|line|
                          {
                              output.append(&mut line.to_vec());
                              output.push(b'\n');
                          }
                         );
                output.append(&mut seq[offset..].to_vec());
                output.push(b'\n');
        }
        println!("{}\n\n{}", String::from_utf8_lossy(fasta), String::from_utf8_lossy(&output));
        assert_eq!(fasta.to_vec(), output);
    }

    #[test]
    fn fastq() {
        let mut x = BufReader::new(Cursor::new(
            "@a\nAGTC\n+\n'&'*+\n@b\nTAGCTTTT\n+\n'&'*+'&'*+\n@c\nGCTA\n+\n'&'*+",
        ));
        let mut record = FastQRecord::default();
        let _ = record.read(&mut x);
        assert_eq!("@a", record.name());
        assert_eq!(b"AGTC".to_vec(), record.seq());
        assert_eq!(&b"AGTC".to_vec(), record.seq_raw());
        let _ = record.read(&mut x);
        assert_eq!("@b", record.name());
        assert_eq!(b"TAGCTTTT".to_vec(), record.seq());
        assert_eq!(&b"TAGCTTTT".to_vec(), record.seq_raw());
        let _ = record.read(&mut x);
        assert_eq!("@c", record.name());
        assert_eq!(b"GCTA".to_vec(), record.seq());
        assert_eq!(&b"GCTA".to_vec(), record.seq_raw());
    }
}
