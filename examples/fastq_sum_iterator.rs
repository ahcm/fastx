use core::iter::Iterator;
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
        let mut len = 0;
        for fastx_record in FastX::fastq_iter(fastx_reader)
        {
            len += fastx_record?.seq_len();
        }
        println!("Total sequence length: {}", len)
    }
    Ok(())
}
