# FastX

FastX implements low overhead readers for Fasta and FastQ.

```
let mut fastx_reader = FastX::reader_from_path(Path::new(&filename))?;
let mut fastx_record = FastX::from_reader(&mut fastx_reader)?;

while let Ok(_some @ 1..=usize::MAX) = fastx_record.read(&mut fastx_reader)
{
  let (id, desc) = fastx_record.name().split_once(" ").unwrap_or((fastx_record.name(), ""));
  println!("{}\t{}\t{}", id, fastx_record.seq_len(), desc)
}


