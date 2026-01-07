// Example: Read a specific chromosome from a remote FASTA file
//
// This example demonstrates fetching chromosome 8 from the Bos taurus (cattle)
// genome hosted on Ensembl's FTP server using indexed random access.
//
// Run with: cargo run --example fasta_fetch_remote_indexed --features url

use std::error::Error;

fn main() -> Result<(), Box<dyn Error>>
{
    // URLs for the data and indexes
    let data_url = "https://ftp.ensembl.org/pub/release-115/fasta/bos_taurus/dna_index/Bos_taurus.ARS-UCD2.0.dna.toplevel.fa.gz";
    let fai_url = "https://ftp.ensembl.org/pub/release-115/fasta/bos_taurus/dna_index/Bos_taurus.ARS-UCD2.0.dna.toplevel.fa.gz.fai";
    let gzi_url = "https://ftp.ensembl.org/pub/release-115/fasta/bos_taurus/dna_index/Bos_taurus.ARS-UCD2.0.dna.toplevel.fa.gz.gzi";

    // Create the remote indexed reader
    let mut reader = fastx::indexed::IndexedFastXReader::from_url(data_url, fai_url, gzi_url)?;

    // Fetch just the first 1000 bases of chromosome 8
    println!("Fetching first 1000 bases of chromosome 8...");
    let seq = reader.fetch_range("8", 0, 1000)?;

    println!("Fetched {} bases", seq.len());
    let first_100 = std::str::from_utf8(&seq[..seq.len().min(100)])
        .unwrap_or("<invalid UTF-8>");
    println!("First 100 bases: {}", first_100);

    Ok(())
}
