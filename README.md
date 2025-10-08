# FastX

FastX implements low overhead readers for Fasta and FastQ.
Version 0.3.0 added .gz support.
Version 0.5.0 added iterator support.

```
let mut fastx_reader = FastX::reader_from_path(Path::new(&filename))?;
let mut fastx_record = FastX::from_reader(&mut fastx_reader)?;

while let Ok(_some @ 1..=usize::MAX) = fastx_record.read(&mut fastx_reader)
{
  println!("{}\t{}", fastx_record.id(), fastx_record.seq_len())
}

```

Or with iterator:

```
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
```
