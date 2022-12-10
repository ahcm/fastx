# FastX

FastX implements low overhead readers for Fasta and FastQ.

```
let mut fastx_reader = FastX::reader_from_path(Path::new(&filename))?;
let mut fastx_record = FastX::from_reader(&mut fastx_reader)?;

while let Ok(_some @ 1..=usize::MAX) = fastx_record.read(&mut fastx_reader)
{
  println!("{}\t{}", fastx_record.id(), fastx_record.seq_len())
}


