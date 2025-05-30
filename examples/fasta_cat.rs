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
        let mut fastx_reader = FastX::reader_from_path(Path::new(&filename))?;
        let mut fastx_record = FastX::FastARecord::default();

        while let Ok(_some @ 1..=usize::MAX) = fastx_record.read(&mut fastx_reader)
        {
            println!(
                ">{} {}\n{}",
                fastx_record.id(),
                fastx_record.desc(),
                String::from_utf8_lossy(&fastx_record.seq())
            );
        }
    }
    Ok(())
}
