use fastx::FastX::open;
use fastx::FastX::FastXRead;
use fastx::FastX::FastXRecord;
use std::env::args;
use std::fs::File;
use std::io;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Read;
use core::iter::Iterator;

fn main() -> io::Result<()>
{
    for filename in args().skip(1) {
        println!("{}", filename);
        let file = File::open(filename).unwrap();
        let mut reader = BufReader::new(file);
        //let mut record = FastXRecord::default();
        let mut record = open(&mut reader)?;
        while let Ok(_some @ 1..=usize::MAX) = record.read(&mut reader)
        {
            let (id, title) = record.name().split_once(" ").unwrap_or((record.name(), ""));
            println!("{}\t{}\t{}", id,
                     record.seq().len(),
                     title)
            //println!("{}\t{}", record.name(), record.seq().len())
        }
    }
    Ok(())
}
