// Example: Read a specific chromosome from a remote FASTA file
//
// This example demonstrates fetching chromosome 8 from the Bos taurus (cattle)
// genome hosted on Ensembl's FTP server using indexed random access.
//
// Run with: cargo run --example remote_chromosome --features url

use std::error::Error;

use fastx::FastX::FastXRead;

fn main() -> Result<(), Box<dyn Error>>
{
    // URLs for the data and indexes
    let data_url = "https://ftp.ensembl.org/pub/release-115/fasta/bos_taurus/dna_index/Bos_taurus.ARS-UCD2.0.dna.toplevel.fa.gz";
    let fai_url = "https://ftp.ensembl.org/pub/release-115/fasta/bos_taurus/dna_index/Bos_taurus.ARS-UCD2.0.dna.toplevel.fa.gz.fai";
    let gzi_url = "https://ftp.ensembl.org/pub/release-115/fasta/bos_taurus/dna_index/Bos_taurus.ARS-UCD2.0.dna.toplevel.fa.gz.gzi";

    // Create the remote indexed reader
    let mut reader = fastx::indexed::IndexedFastXReader::from_url(data_url, fai_url, gzi_url)?;

    // Fetch chromosome 8
    println!("Fetching chromosome 8...");
    let record = reader.fetch("8")?;

    // Get the sequence ID and length
    println!("Sequence ID: {}", record.id());
    println!("Description: {}", record.desc());
    println!("Sequence length: {} bp", record.seq_len());

    Ok(())
}
