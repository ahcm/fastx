// Example: Read a specific chromosome from a FASTA file
//
// Run with: cargo run --example fasta_count_indexed_chromosome seq_id fasta_path

use std::env::args;
use std::path::Path;

use fastx::FastX::FastXRead;

fn main()
{
    let id = args().nth(2).unwrap();
    for filename in args().skip(2)
    {
        println!("{}", filename);

        let data_url = Path::new(&filename);

        let mut reader = fastx::indexed::IndexedFastXReader::from_path(data_url)
            .expect("failed to build Reader from path");

        let record = reader.fetch(&id).expect("record not found");

        println!("{}\t{}", record.id(), record.seq_len());
    }
}
