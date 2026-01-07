# FastX

FastX implements low overhead readers for Fasta and FastQ.
 - Version 0.3.0 added .gz support.
 - Version 0.5.0 added iterator support.
 - Version 0.5.2 added compressed fasta (fasta.gz).
 - Version 0.6.0 added support for URLs and indexed (.fai and .gzi) fasta files.

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

## Features

FastX supports different compression backends through Cargo features. Choose the backend that best fits your needs:

### Default: Pure Rust Backend

By default, FastX uses the `rust-backend` feature, which provides a pure Rust implementation (miniz_oxide) for gzip decompression:

```toml
[dependencies]
fastx = "0.5"
```

This is the safest option with no C dependencies, making it ideal for cross-compilation and environments where you want to avoid native code.

### Alternative Backends

For better performance, you can select different compression backends:

**System zlib** - Uses the zlib library installed on your system (typically fastest on systems with optimized zlib):
```toml
[dependencies]
fastx = { version = "0.5", default-features = false, features = ["zlib"] }
```

**zlib-ng (compatibility mode)** - Modern, fast implementation that's API-compatible with zlib:
```toml
[dependencies]
fastx = { version = "0.5", default-features = false, features = ["zlib-ng-compat"] }
```

**zlib-ng (native API)** - Fastest option using zlib-ng's native API:
```toml
[dependencies]
fastx = { version = "0.5", default-features = false, features = ["zlib-ng"] }
```

### Performance Considerations

- **rust-backend**: Safe, portable, no build dependencies. Moderate performance.
- **zlib**: Good performance if your system has an optimized zlib (e.g., Intel's optimized version).
- **zlib-ng-compat**: Better performance than standard zlib on most systems.
- **zlib-ng**: Best performance, but requires C compiler at build time.

## Documentation
Some people and LLMs might prefer reading documentation instead of examples or code.
So it was generated with claude.
```
● Now I understand the codebase. This is a Rust bioinformatics library for reading FASTA/FASTQ files, and it currently has minimal documentation (on
ly 1 doc comment in the entire codebase). Let me add comprehensive documentation following Rust best practices.
```
