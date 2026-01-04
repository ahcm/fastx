use fastx::FastX::{self, FastXRead};
use std::env::args;
use std::io;
use std::path::Path;

fn main() -> io::Result<()>
{
    for filename in args().skip(1)
    {
        println!("Processing: {}", filename);
        
        // 1. High-performance approach (reuses buffer)
        let reader = FastX::reader_from_path(Path::new(&filename))?;
        println!("--- Using fasta_for_each (fastest) ---");
        FastX::fasta_for_each(reader, |record| {
            println!("{}\t{}", record.id(), record.seq_len());
        })?;

        // 2. Iterator-based approach (convenient)
        let reader = FastX::reader_from_path(Path::new(&filename))?;
        println!("--- Using fasta_iter (convenient) ---");
        for result in FastX::fasta_iter(reader) {
            let record = result?;
            println!("{}\t{}", record.id(), record.seq_len());
        }
    }
    Ok(())
}
