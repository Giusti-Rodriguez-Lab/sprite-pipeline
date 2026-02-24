use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;

// ---------------------------------------------------------------------------
// FASTQ record
// ---------------------------------------------------------------------------

/// A single FASTQ record (4 lines).
///
/// All fields are raw bytes — no UTF-8 assumption.  The header includes the
/// leading `@`.  Newlines are stripped when reading and re-added when writing.
#[derive(Debug, Clone)]
pub struct Record {
    pub header: Vec<u8>,
    pub seq:    Vec<u8>,
    pub plus:   Vec<u8>,
    pub qual:   Vec<u8>,
}

impl Record {
    /// The read name: bytes after `@` up to (but not including) the first
    /// space.  This matches `fq.getName().split(" ")[0]` in the Java version.
    pub fn name_bytes(&self) -> &[u8] {
        let h = self.header.strip_prefix(b"@").unwrap_or(&self.header);
        h.split(|&b| b == b' ')
            .next()
            .unwrap_or(h)
    }

    /// Return a new record whose header is `@<name>::<label>`, where `name`
    /// is the original read name (no `@`, no trailing comment) and `label`
    /// is the concatenated barcode string built by the matcher.
    ///
    /// Mirrors Java's `appendBarcodesToName`:
    /// `fq.changeName(fq.getName().split(" ")[0] + "::" + sb.toString())`
    pub fn with_barcode_label(&self, label: &[u8]) -> Self {
        let name = self.name_bytes();
        let mut header: Vec<u8> = Vec::with_capacity(1 + name.len() + 2 + label.len());
        header.push(b'@');
        header.extend_from_slice(name);
        header.extend_from_slice(b"::");
        header.extend_from_slice(label);

        Record {
            header,
            seq:  self.seq.clone(),
            plus: self.plus.clone(),
            qual: self.qual.clone(),
        }
    }

    /// Serialize the record as 4 newline-terminated lines.
    pub fn write_to<W: Write>(&self, w: &mut W) -> std::io::Result<()> {
        w.write_all(&self.header)?;
        w.write_all(b"\n")?;
        w.write_all(&self.seq)?;
        w.write_all(b"\n")?;
        w.write_all(&self.plus)?;
        w.write_all(b"\n")?;
        w.write_all(&self.qual)?;
        w.write_all(b"\n")?;
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Streaming FASTQ reader
// ---------------------------------------------------------------------------

pub struct FastqReader<R: BufRead> {
    reader: R,
    buf:    Vec<u8>,
}

impl<R: BufRead> FastqReader<R> {
    pub fn new(reader: R) -> Self {
        Self { reader, buf: Vec::with_capacity(1024) }
    }

    /// Read the next record.  Returns `Ok(None)` at end-of-file.
    pub fn next_record(&mut self) -> std::io::Result<Option<Record>> {
        // --- header ---
        self.buf.clear();
        let n = self.reader.read_until(b'\n', &mut self.buf)?;
        if n == 0 {
            return Ok(None); // clean EOF
        }
        let header = trim_line_ending(&self.buf);
        if header.is_empty() {
            return Ok(None);
        }

        // --- sequence ---
        self.buf.clear();
        self.reader.read_until(b'\n', &mut self.buf)?;
        let seq = trim_line_ending(&self.buf);

        // --- '+' line ---
        self.buf.clear();
        self.reader.read_until(b'\n', &mut self.buf)?;
        let plus = trim_line_ending(&self.buf);

        // --- quality ---
        self.buf.clear();
        self.reader.read_until(b'\n', &mut self.buf)?;
        let qual = trim_line_ending(&self.buf);

        Ok(Some(Record { header, seq, plus, qual }))
    }
}

fn trim_line_ending(buf: &[u8]) -> Vec<u8> {
    let end = buf.len();
    if end > 0 && buf[end - 1] == b'\n' {
        if end > 1 && buf[end - 2] == b'\r' {
            buf[..end - 2].to_vec()
        } else {
            buf[..end - 1].to_vec()
        }
    } else {
        buf.to_vec()
    }
}

// ---------------------------------------------------------------------------
// Convenience constructors
// ---------------------------------------------------------------------------

/// Open a gzipped FASTQ file for reading.
pub fn open_gz(path: &Path) -> std::io::Result<FastqReader<BufReader<GzDecoder<File>>>> {
    let file    = File::open(path)?;
    let decoder = GzDecoder::new(file);
    Ok(FastqReader::new(BufReader::with_capacity(1 << 20, decoder)))
}

/// Create a gzipped FASTQ output file.
pub fn create_gz(path: &Path) -> std::io::Result<BufWriter<GzEncoder<File>>> {
    let file    = File::create(path)?;
    let encoder = GzEncoder::new(file, Compression::default());
    Ok(BufWriter::with_capacity(1 << 20, encoder))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn make_reader(s: &str) -> FastqReader<Cursor<&[u8]>> {
        FastqReader::new(Cursor::new(s.as_bytes()))
    }

    #[test]
    fn test_reads_single_record() {
        let fq = "@read1\nACGT\n+\nIIII\n";
        let mut reader = make_reader(fq);
        let rec = reader.next_record().unwrap().unwrap();
        assert_eq!(rec.header, b"@read1");
        assert_eq!(rec.seq,    b"ACGT");
        assert_eq!(rec.plus,   b"+");
        assert_eq!(rec.qual,   b"IIII");
    }

    #[test]
    fn test_reads_two_records() {
        let fq = "@r1\nAAAA\n+\nIIII\n@r2\nCCCC\n+\nHHHH\n";
        let mut reader = make_reader(fq);
        let r1 = reader.next_record().unwrap().unwrap();
        let r2 = reader.next_record().unwrap().unwrap();
        assert_eq!(r1.header, b"@r1");
        assert_eq!(r2.header, b"@r2");
    }

    #[test]
    fn test_eof_returns_none() {
        let fq = "@r1\nAAAA\n+\nIIII\n";
        let mut reader = make_reader(fq);
        reader.next_record().unwrap().unwrap();
        let eof = reader.next_record().unwrap();
        assert!(eof.is_none());
    }

    #[test]
    fn test_name_bytes_strips_at_and_comment() {
        let rec = Record {
            header: b"@SRR123.1 comment here".to_vec(),
            seq:    vec![],
            plus:   vec![],
            qual:   vec![],
        };
        assert_eq!(rec.name_bytes(), b"SRR123.1");
    }

    #[test]
    fn test_name_bytes_no_comment() {
        let rec = Record {
            header: b"@SRR123.1".to_vec(),
            seq:    vec![],
            plus:   vec![],
            qual:   vec![],
        };
        assert_eq!(rec.name_bytes(), b"SRR123.1");
    }

    #[test]
    fn test_with_barcode_label() {
        let rec = Record {
            header: b"@SRR123.1 comment".to_vec(),
            seq:    b"ACGT".to_vec(),
            plus:   b"+".to_vec(),
            qual:   b"IIII".to_vec(),
        };
        let labelled = rec.with_barcode_label(b"[DPM_A][Y_1]");
        assert_eq!(labelled.header, b"@SRR123.1::[DPM_A][Y_1]");
        // Sequence and quality untouched
        assert_eq!(labelled.seq,  b"ACGT");
        assert_eq!(labelled.qual, b"IIII");
    }

    #[test]
    fn test_write_to_roundtrip() {
        let fq = "@read1\nACGT\n+\nIIII\n";
        let mut reader = make_reader(fq);
        let rec = reader.next_record().unwrap().unwrap();

        let mut out = Vec::new();
        rec.write_to(&mut out).unwrap();
        assert_eq!(out, fq.as_bytes());
    }

    #[test]
    fn test_crlf_line_endings() {
        let fq = "@r1\r\nACGT\r\n+\r\nIIII\r\n";
        let mut reader = make_reader(fq);
        let rec = reader.next_record().unwrap().unwrap();
        assert_eq!(rec.header, b"@r1");
        assert_eq!(rec.seq,    b"ACGT");
    }
}
