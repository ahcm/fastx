use criterion::{black_box, criterion_group, criterion_main, Criterion};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use rand::Rng;
use fastx::FastX::FastXRead;

fn generate_fasta(path: &Path, size_mb: usize) {
    let mut file = BufWriter::new(File::create(path).unwrap());
    let mut rng = rand::thread_rng();
    let bases = b"ACGT";
    let line_len = 80;
    
    let mut written = 0;
    let target = size_mb * 1024 * 1024;
    let mut i = 0;
    
    while written < target {
        writeln!(file, ">seq{}", i).unwrap();
        written += 10; // Approx header len
        
        let seq_len = rng.gen_range(100..1000);
        for j in 0..seq_len {
            file.write_all(&[bases[rng.gen_range(0..4)]]).unwrap();
            if (j + 1) % line_len == 0 {
                file.write_all(b"\n").unwrap();
            }
        }
        file.write_all(b"\n").unwrap();
        written += seq_len;
        i += 1;
    }
}

fn bench_fastx(c: &mut Criterion) {
    let file_path = Path::new("bench_data.fasta");
    if !file_path.exists() {
        generate_fasta(file_path, 10);
    }

    let mut group = c.benchmark_group("parsing");

    group.bench_function("fastx iterator", |b| {
        b.iter(|| {
            let reader = std::io::BufReader::new(File::open(file_path).unwrap());
            let mut count = 0;
            let mut bases = 0;
            for result in fastx::FastX::fasta_iter(reader) {
                let record = result.unwrap();
                count += 1;
                bases += record.seq_len();
                black_box(record.id());
            }
            black_box((count, bases));
        })
    });

    group.bench_function("fastx manual reuse", |b| {
        b.iter(|| {
            let mut reader = std::io::BufReader::new(File::open(file_path).unwrap());
            let mut record = fastx::FastX::FastARecord::default();
            let mut count = 0;
            let mut bases = 0;
            while let Ok(n) = record.read(&mut reader) {
                if n == 0 { break; }
                count += 1;
                bases += record.seq_len();
                black_box(record.id());
            }
            black_box((count, bases));
        })
    });

    group.bench_function("fastx for_each", |b| {
        b.iter(|| {
            let mut reader = std::io::BufReader::new(File::open(file_path).unwrap());
            let mut record = fastx::FastX::FastARecord::default();
            let mut count = 0;
            let mut bases = 0;
            let mut func = |rec: &fastx::FastX::FastARecord| {
                count += 1;
                bases += rec.seq_len();
                black_box(rec.id());
            };
            
            while let Ok(n) = record.read(&mut reader) {
                if n == 0 { break; }
                func(&record);
            }
            black_box((count, bases));
        })
    });

    group.bench_function("needletail", |b| {
        b.iter(|| {
            let mut reader = needletail::parse_fastx_file(file_path).unwrap();
            let mut count = 0;
            let mut bases = 0;
            while let Some(record) = reader.next() {
                let record = record.unwrap();
                count += 1;
                bases += record.num_bases();
                black_box(record.id());
            }
            black_box((count, bases));
        })
    });
    
    group.finish();
}

fn bench_needletail(c: &mut Criterion) {
    // Merged into bench_fastx for better grouping
}

criterion_group!(benches, bench_fastx);
criterion_main!(benches);
