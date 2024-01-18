use core::iter::Iterator;
use fastx::FastX;
use std::env::args;
use std::io;
use std::path::Path;

fn main() -> io::Result<()>
{
    for filename in args().skip(1)
    {
        println!("{}", filename);
        let mut fastx_reader = FastX::reader_from_path(Path::new(&filename))?;
        let mut fastx_record = FastX::from_reader(&mut fastx_reader)?;

        /*
        // just for fun read the first line and seek back
        use std::io::BufRead;
        let mut line = String::new();
        let offset = fastx_reader.read_line(&mut line)?;
        println!("{}", line);
        fastx_reader.seek_relative(- (offset as i64))?;
        */

        // back to serious processing
        while let Ok(_some @ 1..=usize::MAX) = fastx_record.read(&mut fastx_reader)
        {
            println!("{}\t{}", fastx_record.id(), fastx_record.seq_len())
        }
    }
    Ok(())
}
