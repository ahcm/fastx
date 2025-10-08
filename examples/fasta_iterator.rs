use fastx::FastX::{self, FastXRead};
use std::env::args;
use std::io;
use std::path::Path;

fn main() -> io::Result<()>
{
    for filename in args().skip(1)
    {
        println!("{}", filename);
        let fastx_reader = FastX::reader_from_path(Path::new(&filename))?;

        // Using the new iterator approach
        for result in FastX::fasta_iter(fastx_reader) {
            match result {
                Ok(record) => println!("{}\t{}", record.id(), record.seq_len()),
                Err(e) => eprintln!("Error reading record: {}", e),
            }
        }
    }
    Ok(())
}
