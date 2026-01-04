# FastX

FastX implements low overhead readers for Fasta and FastQ.
 - Version 0.3.0 added .gz support.
 - Version 0.5.0 added iterator support.
 - Version 0.5.2 added compressed fasta (fasta.gz).
 - Version 0.6.0 added high-performance `for_each` iteration, improved FASTQ parsing correctness, and exposed public API for indexed random access (.fai and .gzi).

## Examples

### High-Performance Reading (Recommended)

For maximum performance, use the `for_each` methods which reuse the internal buffer for each record, avoiding allocations.

```rust
use fastx::FastX::{self, FastXRead};
use std::path::Path;

fn main() -> std::io::Result<()> {
    let path = Path::new("sequences.fastq.gz");
    let reader = FastX::reader_from_path(path)?;

    // fastx_for_each automatically detects FASTA vs FASTQ
    FastX::fastx_for_each(reader,
        |record| {
            println!("FASTA: {}\tLen: {}", record.id(), record.seq_len());
        },
        |record| {
            println!("FASTQ: {}\tLen: {}", record.id(), record.seq_len());
        }
    )?;
    Ok(())
}
```

### Iterator-Based Reading (Convenient)

The iterator API is convenient for standard Rust loops but allocates a new buffer for each record.

```rust
use fastx::FastX::{self, FastXRead};
use std::path::Path;

fn main() -> std::io::Result<()> {
    let path = Path::new("sequences.fasta");
    let reader = FastX::reader_from_path(path)?;

    for result in FastX::fasta_iter(reader) {
        let record = result?;
        println!("{}\t{}", record.id(), record.seq_len());
    }
    Ok(())
}
```

### Random Access with Indexed Files

FastX supports random access to BGZF-compressed FASTA files using `.fai` and `.gzi` indexes.

```rust
use fastx::indexed::IndexedFastXReader;
use fastx::FastX::FastXRead;
use std::path::Path;

fn main() -> std::io::Result<()> {
    // Requires data.fasta.gz, data.fasta.gz.fai, and data.fasta.gz.gzi
    let mut reader = IndexedFastXReader::from_path(Path::new("data.fasta.gz"))?;

    // Fetch a specific sequence by ID
    let record = reader.fetch("chr1")?;
    println!("{}: {} bp", record.id(), record.seq_len());

    // Fetch a specific region (0-based start, exclusive end)
    let region = reader.fetch_range("chr1", 1000, 2000)?;
    println!("Region length: {}", region.len());
    
    Ok(())
}
```

## Features

FastX supports different compression backends through Cargo features. Choose the backend that best fits your needs:

### Default: Pure Rust Backend

By default, FastX uses the `rust-backend` feature, which provides a pure Rust implementation (miniz_oxide) for gzip decompression:

```toml
[dependencies]
fastx = "0.6"
```

This is the safest option with no C dependencies, making it ideal for cross-compilation and environments where you want to avoid native code.

### Alternative Backends

For better performance, you can select different compression backends:

**System zlib** - Uses the zlib library installed on your system (typically fastest on systems with optimized zlib):
```toml
[dependencies]
fastx = { version = "0.6", default-features = false, features = ["zlib"] }
```

**zlib-ng (compatibility mode)** - Modern, fast implementation that's API-compatible with zlib:
```toml
[dependencies]
fastx = { version = "0.6", default-features = false, features = ["zlib-ng-compat"] }
```

**zlib-ng (native API)** - Fastest option using zlib-ng's native API:
```toml
[dependencies]
fastx = { version = "0.6", default-features = false, features = ["zlib-ng"] }
```

### URL Support

Enable the `url` feature (enabled by default) to read indexed files directly from HTTP/HTTPS URLs:

```rust
let mut reader = IndexedFastXReader::from_url(
    "https://example.com/data.fasta.gz",
    "https://example.com/data.fasta.gz.fai",
    "https://example.com/data.fasta.gz.gzi"
)?;
```

### Performance Considerations

- **rust-backend**: Safe, portable, no build dependencies. Moderate performance.
- **zlib**: Good performance if your system has an optimized zlib (e.g., Intel's optimized version).
- **zlib-ng-compat**: Better performance than standard zlib on most systems.
- **zlib-ng**: Best performance, but requires C compiler at build time.