use fastx::FastX::{self, FastQRead};
use fastx::FastX::FastXRead;
use std::env::args;
use std::io;
use std::path::Path;
use core::iter::Iterator;

fn main() -> io::Result<()>
{
    for filename in args().skip(1)
    {
        println!("{}", filename);
        let mut fastx_reader = FastX::reader_from_path(Path::new(&filename))?;
        let mut fastq_record = FastX::FastQRecord::default();

        while let Ok(_some @ 1..=usize::MAX) = fastq_record.read(&mut fastx_reader)
        {
            println!("@{} {}\n{}\n+{}\n{}",
                     fastq_record.id(),
                     fastq_record.desc(),
                     String::from_utf8_lossy(&fastq_record.seq()),
                     fastq_record.comment(),
                     String::from_utf8_lossy(&fastq_record.qual()));
        }
    }
    Ok(())
}
