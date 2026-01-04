use fastx::FastX::{self, FastXRead};
use std::env::args;
use std::io;
use std::path::Path;

fn main() -> io::Result<()>
{
    for filename in args().skip(1)
    {
        println!("File: {}", filename);
        let reader = FastX::reader_from_path(Path::new(&filename))?;

        // fastx_for_each automatically detects if it is FASTA or FASTQ
        // and uses the appropriate high-performance loop.
        FastX::fastx_for_each(reader,
            |record| {
                // This closure is called for FASTA records
                println!("FASTA ID: {}\tLength: {}", record.id(), record.seq_len());
            },
            |record| {
                // This closure is called for FASTQ records
                println!("FASTQ ID: {}\tLength: {}", record.id(), record.seq_len());
            }
        )?;
    }
    Ok(())
}
