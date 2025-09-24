
"""
Figure out the input file to use and tell the user if they have an invalid file.
"""

# ------------------------- SB_comments (A): background -------------------------
# # Looking at the information at: https://ablab.github.io/spades/input.html, and searching the internet:
#
# (1) SPAdes directly accepts the following file formats for input:
#     ['.fastq', '.fq', '.fastq.gz', '.fq.gz']  (primary read formats)
#     ['.fna', '.fa', '.fasta']                 (for trusted/untrusted contigs, references)
#
# (2) For alignment files: ['.bam', '.sam', '.cram'], 
#     they must first be converted to FASTQ using samtools:
#         e.g., samtools fastq input.bam -1 reads_1.fastq -2 reads_2.fastq -s singles.fastq
#
# (3) For SRA archive files: ['.sra'], 
#     they must first be extracted to FASTQ using the SRA Toolkit (fastq-dump or fasterq-dump):
#         e.g., fasterq-dump input.sra -O ./output_folder 
#         (produces files like SRR1234567_1.fastq and SRR1234567_2.fastq)
#
# (4) For PacBio legacy formats: ['.h5', '.bax.h5'], 
#     they must first be converted to PacBio BAM, then to FASTQ:
#         Step 1: bax2bam movie.bax.h5 -o movie.subreads.bam
#         Step 2: bam2fastq movie.subreads.bam -o movie_reads
#
# (5) For Oxford Nanopore raw signal files: ['.fast5'], 
#     they must be basecalled to FASTQ using Guppy or Dorado:
#         e.g., guppy_basecaller -i input_fast5/ -s output_folder/ --flowcell FLO-MIN106 --kit SQK-LSK109
#         e.g., dorado basecaller dna_r9.4.1 -r input_fast5/ > output.fastq
#     (Both tools can also output BAM if requested: Guppy --bam_out, Dorado --emit-bam)

##STRATEGY
'''
I deployed the following Strategy

1) First accept a path from the command line (either a folder or a single file). I did this so
   people did not need to edit code to point at their data.

2) Then, list all files in a clean, numbered catalog. For each file, print its extension, a short preview
   (first lines for text, hex signature for binaries), and a quick content-based classification
   (FASTQ, FASTA, SAM, BAM, CRAM, HDF5/FAST5, SRA, or unknown). I did this so users could
   quickly see “what they actually had” before assembly.

3) Ask the user to select items by number and then choose an action for each (use / convert / skip).
   I used “convert” to print precise hints for turning alignments/containers into FASTQ, because SPAdes
   ultimately wants reads.

4) For FASTQ files that were marked “use,” pair R1 with R2 and ran a small sanity check:
   by a handful of normalized read IDs to estimate how much the two files overlapped. I did this an an extra feature
   so that users can catch wrong pairings before SPAdes.

5) Print a SPAdes command template and a compact options cheat-sheet.

6) I directed the script to write a CSV report (catalog + R1/R2 pairing sanity) for record-keeping.

I have explained in line, why I chose each construct and which import or built-in packages I used.
'''

import sys
    # access to Python runtime info and command-line arguments.
    # Read command-line arguments (e.g., the user’s folder path).
    # Exit with status codes (sys.exit(1) for errors).
    # Know which Python we’re running on (sys.version).

import os
    # operating system interface (files, environment variables, permissions).
    # Check file readability (os.access(path, os.R_OK)).
    # Work with directories (os.listdir, os.makedirs).
    # Read environment variables (os.environ).

import re
    # regular expressions: advanced pattern matching in text.
    # Detect R1/R2 or “forward/reverse” in filenames.
    # Clean or normalize read headers (remove /1, /2 suffix).
    # Flexible string searching/replacing.

import gzip
    # open compressed .gz files (common for FASTQ/FASTA).
    # Treat gzipped files like normal text or binary streams.
    # Example: gzip.open("reads.fastq.gz", "rt") for reading as text.

import csv
    # read and write CSV (comma-separated values) files.
    # Save reports and tables in spreadsheet-friendly format.
    # Example: writer.writerow([...]) to save one row in a CSV.

import binascii
    # convert binary data to hex string.
    # Preview "magic bytes" of binary files (e.g., BAM, CRAM, HDF5).
    # Useful for format detection when file is not plain text.

import argparse
    # command-line argument parsing.
    # Define user flags and options (e.g., --show-others).
    # Automatically generates --help message.
    # Cleaner and safer than manually handling sys.argv.

from pathlib import Path
    # object-oriented file path handling.
    # Cross-platform and safer than os.path.
    # Examples:
    #   Path("folder") / "reads.fastq.gz" → folder/reads.fastq.gz
    #   path.exists(), path.is_file(), path.suffixes, path.iterdir()


# I set up a few small helpers and constants to keep later code short and readable.
# I used these so users could follow the main flow without diving into low-level details.

bold_font, reset_font = "\033[1m", "\033[0m"    # I used ANSI escape codes to make headings stand out in the terminal.
def make_bold(text):
    # I wrapped strings in the bold codes above; this was cosmetic but improved readability of the catalog.
    return f"{bold_font}{text}{reset_font}"

# I defined the sets below because SPAdes either used these directly (FASTQ/FASTA)
# or needed them converted first (alignments, SRA, HDF5/FAST5).
# I used sets so I could test membership clearly and quickly.

fastq_files      = {".fastq", ".fq", ".fastq.gz", ".fq.gz"}
fasta_files      = {".fasta", ".fa", ".fna", ".fasta.gz", ".fa.gz", ".fna.gz", ".faa", ".faa.gz"}
alignment_files  = {".bam", ".sam", ".cram"}
pacbio_files     = {".bax.h5", ".h5"}
nanopore_files   = {".fast5"}
sra_files        = {".sra"}

# I combined everything into a single whitelist of “needed” extensions.
# I used this to keep the file list focused by default.
needed_extensions = set().union(
    fastq_files,
    fasta_files,
    alignment_files,
    pacbio_files,
    nanopore_files,
    sra_files,
)



def all_suffixes_as_one(path_obj: Path) -> str:
    '''
    I wrote this to return a lowercase “extension chain” like ".fastq.gz" (or ".dir" for folders).
    I did this because double-suffix formats were common in bioinformatics and using only the last suffix was misleading.
    I used Path.is_dir and Path.suffixes from pathlib, and I used built-in string join/lower.
    '''
    return ".dir" if path_obj.is_dir() else "".join(path_obj.suffixes).lower()


def open_as_text(path_obj: Path):
    '''
    I wrote this helper to open files as text whether they were gzipped or not.
    I did this so users did not need to gunzip reads manually.
    I used gzip.open in text mode for .gz, and the built-in open for everything else.
    '''
    return gzip.open(path_obj, "rt", encoding="utf-8", errors="replace") \
        if all_suffixes_as_one(path_obj).endswith(".gz") \
        else open(path_obj, "rt", encoding="utf-8", errors="replace")


def open_as_binary(path_obj: Path):
    '''
    I wrote this to open files as binary streams when I only needed the first bytes.
    I did this so I could print concise hex “signatures” for binary formats like BAM/CRAM/HDF5/SRA.
    I used gzip.open in binary mode for .gz and built-in open for plain files.
    '''
    return gzip.open(path_obj, "rb") if all_suffixes_as_one(path_obj).endswith(".gz") else open(path_obj, "rb")


def show_preview(path_obj: Path, n_lines=6, n_bytes=48) -> str:
    '''
    I wrote this to show a quick preview for each file. I used the first few text lines for FASTQ/FASTA/SAM
    and a short hex signature for binaries. I did this to give users a visual check of “does this look right?”
    without overwhelming them.

    I used open_as_text/open_as_binary (my helpers), binascii.hexlify for hex, and small try/except blocks for clarity.
    '''
    ext = all_suffixes_as_one(path_obj)
    if ext in fastq_files | fasta_files | {".sam"}:
        try:
            lines = []
            with open_as_text(path_obj) as fh:
                for _ in range(n_lines):
                    line = fh.readline()
                    if not line:
                        break
                    line = line.rstrip("\n")
                    # I trimmed very long lines so the preview stayed inside typical terminal widths.
                    lines.append(line[:160] + (" …" if len(line) > 160 else ""))
            return "\n".join(lines) or "[empty file]"
        except Exception as err:
            return f"[text preview failed: {err}]"
    try:
        with open_as_binary(path_obj) as fh:
            head = fh.read(n_bytes)
        hx = binascii.hexlify(head).decode("ascii")
        return " ".join(hx[i:i+2] for i in range(0, len(hx), 2)) or "[empty file]"
    except Exception as err:
        return f"[hex preview failed: {err}]"


def looks_like_fastq(path_obj: Path) -> bool:
    '''
    I wrote this to perform a small, structural FASTQ check. I read a few 4-line records and verified:
    1) the first line started with "@",
    2) the third line started with "+",
    3) the sequence and quality lines were equal length.

    I did this because mislabeled files can be common and this fast check prevented surprises downstream.
    I used my open_as_text helper and built-in string methods.
    '''
    try:
        with open_as_text(path_obj) as fh:
            for _ in range(10):
                h = fh.readline(); s = fh.readline(); plus = fh.readline(); q = fh.readline()
                if not (h and s and plus and q):
                    return False
                if not h.startswith("@") or not plus.startswith("+"):
                    return False
                if len(s.rstrip()) != len(q.rstrip()):
                    return False
        return True
    except Exception:
        return False


def looks_like_fasta(path_obj: Path) -> bool:
    '''
    I wrote this to quickly identify FASTA/FAA files. I scanned only a few hundred lines and
    checked for at least one “>” header and at least some non-empty sequence lines.
    I did this because users often mixed FASTA contigs with reads, and I wanted to label them clearly.

    I used open_as_text and enumerate, which kept the code readable.
    '''
    try:
        headers = seq_lines = 0
        with open_as_text(path_obj) as fh:
            for i, line in enumerate(fh, 1):
                if line.startswith(">"):
                    headers += 1
                elif line.strip():
                    seq_lines += 1
                if i >= 400:
                    break
        return headers > 0 and seq_lines > 0
    except Exception:
        return False


def starts_with_magic(path_obj: Path, sig: bytes) -> bool:
    '''
    I wrote this to detect binary formats by their “magic” bytes:
    - BAM began with BAM\\x01
    - CRAM began with CRAM
    - HDF5/FAST5 had a known signature

    I did this because a short binary comparison was much faster and clearer than full parsing.
    I used open_as_binary and a simple equality check.
    '''
    try:
        with open_as_binary(path_obj) as fh:
            return fh.read(len(sig)) == sig
    except Exception:
        return False


def classify_by_look(path_obj: Path):
    '''
    I wrote this classifier so I could give immediate, practical advice (“use this as reads”
    vs “convert this first”). I moved from simple checks to slightly more specific ones:

    1) I tried looks_like_fastq and looks_like_fasta. If either passed, I marked it VALID for SPAdes input.
    2) I checked BAM/CRAM/HDF5 magic bytes; these were NOT VALID for direct SPAdes reads.
    3) I sniffed SAM by looking for header lines or 11+ tab columns in the first ~100 lines.
    4) I looked for “NCBI.sra” near the head for SRA archives.

    I returned a tuple: (kind, is_valid, short_reason), where kind was one of
    fastq, fasta, bam, cram, sam, hdf5, sra, unknown.
    '''
    if looks_like_fastq(path_obj):
        return ("fastq", True,  "FASTQ 4-line records looked correct")
    if looks_like_fasta(path_obj):
        return ("fasta", True,  "FASTA headers and sequence lines were present")
    if starts_with_magic(path_obj, b"BAM\x01"):
        return ("bam",   False, "BAM signature was present; convert to FASTQ first")
    if starts_with_magic(path_obj, b"CRAM"):
        return ("cram",  False, "CRAM signature was present; convert to FASTQ with a reference")
    seems_sam = False
    try:
        with open_as_text(path_obj) as fh:
            for i, line in enumerate(fh, 1):
                if not line.strip():
                    continue
                if line.startswith("@") or len(line.rstrip("\n").split("\t")) >= 11:
                    seems_sam = True
                    break
                if i >= 100:
                    break
    except Exception:
        pass
    if seems_sam:
        return ("sam", False, "Text looked like SAM; convert to FASTQ")
    if starts_with_magic(path_obj, b"\x89HDF\r\n\x1a\n"):
        return ("hdf5", False, "HDF5/FAST5 container detected; basecall/convert to FASTQ")
    try:
        with open_as_binary(path_obj) as fh:
            if b"NCBI.sra" in fh.read(512):
                return ("sra", False, "SRA archive marker found; extract FASTQ")
    except Exception:
        pass
    return ("unknown", False, "I could not identify the type with quick checks")


# I added the following helpers to guess R1/R2 from both names and headers,
# and to verify that R1 and R2 probably belonged to the same sample.

name_has_r1_pattern = re.compile(r"(?:^|[_\.-])(R?1|read1|forward)(?:[^A-Za-z0-9]|$)", re.I)
name_has_r2_pattern = re.compile(r"(?:^|[_\.-])(R?2|read2|reverse)(?:[^A-Za-z0-9]|$)", re.I)
head_has_r1_hint    = re.compile(r"(?:\s|/)(1)(?::|$)", re.I)
head_has_r2_hint    = re.compile(r"(?:\s|/)(2)(?::|$)", re.I)

def guess_orientation_from_name(path_obj: Path) -> str:
    '''
    I used filename patterns first because they were fast and often accurate.
    I returned "R1", "R2", or "?" if the name did not clearly indicate one side.
    '''
    nm = path_obj.name
    if name_has_r1_pattern.search(nm) and not name_has_r2_pattern.search(nm):
        return "R1"
    if name_has_r2_pattern.search(nm) and not name_has_r1_pattern.search(nm):
        return "R2"
    return "?"


def guess_orientation_from_headers(path_obj: Path, max_reads=60) -> str:
    '''
    I fell back to header tags when names were ambiguous. I sampled a few records and counted
    whether “1” or “2” markers showed up more often. If one side clearly dominated, I reported it.
    '''
    r1 = r2 = 0
    try:
        with open_as_text(path_obj) as fh:
            for _ in range(max_reads):
                h = fh.readline()
                if not h:
                    break
                fh.readline(); fh.readline(); fh.readline()
                if h.startswith("@"):
                    if head_has_r1_hint.search(h):
                        r1 += 1
                    if head_has_r2_hint.search(h):
                        r2 += 1
    except Exception:
        pass
    if r1 > r2 and r1 >= 5:
        return "R1"
    if r2 > r1 and r2 >= 5:
        return "R2"
    return "?"


def decide_orientation(path_obj: Path) -> str:
    '''
    I centralized the orientation logic in this function so later code stayed clean.
    I tried name-based detection first and used headers only if names were unclear.
    '''
    by_name = guess_orientation_from_name(path_obj)
    return by_name if by_name != "?" else guess_orientation_from_headers(path_obj)


def normalize_read_id(header_line: str) -> str:
    '''
    I normalized raw FASTQ headers into comparable IDs. I did this by:
    1) stripping whitespace,
    2) keeping only the first whitespace-separated token,
    3) removing a trailing “/1”, “/2”, “:1” or “:2”.
    This let me compare R1 and R2 read sets to estimate overlap.
    '''
    txt = header_line.strip()
    if " " in txt:
        txt = txt.split(" ", 1)[0]
    txt = re.sub(r"/[12]$", "", txt)
    txt = re.sub(r":[12]$", "", txt)
    return txt


def first_read_ids(path_obj: Path, n=100) -> set:
    '''
    I collected a small set of normalized IDs from each FASTQ file.
    I did this so I could estimate overlap between R1 and R2 by set intersection,
    which was a quick and robust sanity check.
    '''
    ids = set()
    try:
        with open_as_text(path_obj) as fh:
            for _ in range(n):
                h = fh.readline()
                if not h:
                    break
                s = fh.readline(); fh.readline(); fh.readline()
                if not s:
                    break
                if h.startswith("@"):
                    ids.add(normalize_read_id(h))
    except Exception:
        pass
    return ids


def infer_sample_stub(path_obj: Path) -> str:
    '''
    I inferred a “sample stub” by removing .fastq/.fq (with optional .gz) and stripping tokens like
    R1, R2, read1, read2, forward, reverse from the filename. I did this grouping so that
    sample_X_R1.fastq.gz naturally paired with sample_X_R2.fastq.gz.
    '''
    name = path_obj.name
    for suf in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.lower().endswith(suf):
            name = name[:-len(suf)]
    name = re.sub(r"([._-])(R?1|R?2|read1|read2|forward|reverse)([^A-Za-z0-9]|$)", r"\1", name, flags=re.I)
    return name.strip("._-")


def spades_template_line() -> str:
    '''
    I printed a template instead of a full command because the ideal options depended on the dataset
    and the user’s goal. The template reminded them of the overall structure.
    '''
    return "spades.py <options (e.g., -1 <R1> -2 <R2, etc)>  (see below for more details)> -o <output>"


def spades_options_block() -> str:
    '''
    I included a short options list here so people did not need to leave the script to recall
    common SPAdes flags. I kept it small and practical.
    '''
    return """\
[SPAdes options you can append in <options>]

General Options
These flags control the core behavior, performance, and resource usage of the SPAdes assembler.

        --careful

            What it is: This flag enables an additional, computationally intensive mismatch and short-indel correction algorithm called MismatchCorrector after the initial BayesHammer error correction.

            Critical Info: It significantly increases runtime and memory usage but can produce a higher-fidelity assembly with fewer small errors (SNPs and indels). It's highly recommended for isolate genomes where accuracy is paramount (e.g., for SNP analysis) but is often too slow for large, complex metagenomes.

        --threads <N>

            What it is: Sets the number of CPU threads to use for parallel computation.

            Critical Info: The more threads you provide, the faster the assembly will run, up to a point. Performance doesn't always scale linearly. A safe bet is to use the number of available cores on the machine (e.g., --threads 16).

        --memory <GB>

            What it is: Specifies the maximum amount of RAM in Gigabytes that SPAdes is allowed to use.

            Critical Info: This is a crucial parameter to prevent the assembly job from crashing. SPAdes can be memory-hungry, especially with large k-mer sizes and complex datasets. Always set this to a value slightly less than the total available RAM on the system (e.g., --memory 50 on a 64 GB machine).

        -k <k1,k2,...>

            What it is: Specifies the k-mer sizes for the de Bruijn graph. SPAdes uses an iterative approach, starting with the smallest k-mer and proceeding to the largest.

            Critical Info: This is arguably the most important parameter affecting assembly contiguity.

            Small k-mers (e.g., 21, 33) are robust to low coverage but struggle to resolve repeats.

            Large k-mers (e.g., 99, 127) are excellent at resolving repeats but are more sensitive to low coverage and errors.

            By using multiple k-mers, SPAdes gets the best of both worlds. The default (-k 21,33,55,77,99,127) is a good starting point for Illumina data on bacterial genomes. Crucially, the maximum k-mer size must be an odd number and smaller than the read length.

        --cov-cutoff <auto|float>

            What it is: Sets a threshold to remove low-coverage edges from the assembly graph. These edges are often due to sequencing errors.

            Critical Info: The default auto setting works well in most cases. However, for metagenomes or isolate genomes where you suspect the presence of low-coverage plasmids, manually setting a lower cutoff (e.g., --cov-cutoff 2.0) might help recover them. Conversely, a higher cutoff can help clean up noisy data.

        -o <outdir>

            What it is: The directory where all output files will be saved.

            Critical Info: This parameter is mandatory. SPAdes requires that the specified output directory does not already exist. This is a safety feature to prevent you from accidentally overwriting a previous assembly.

    Data-type Presets
    These are essential flags that tune numerous internal parameters to optimize the assembly algorithm for a specific type of data.

        --isolate

            What it is: The default mode, optimized for assembling the genome of a single organism (a pure culture) where coverage is expected to be relatively uniform.

        --metagenome

            What it is: A specialized mode for assembling complex, mixed samples containing many different species with varying abundances (e.g., gut or soil samples).

            Critical Info: This mode uses different graph simplification algorithms and coverage models to handle the high complexity and uneven coverage inherent in metagenomic data. Do not use this for isolate data.

        --rna

            What it is: Invokes the rnaSPAdes pipeline for transcriptome assembly from RNA-Seq data. It is designed to handle isoform complexity.

        --sc

            What it is: Invokes the scSPAdes pipeline, which is specifically designed for single-cell sequencing data. This data is characterized by extremely non-uniform coverage due to whole-genome amplification biases, which this mode is built to handle.

    Reads
    These flags specify the input sequencing data. The files are typically in FASTQ format (.fastq or .fq), often compressed with gzip (.gz).

        -1 <R1> -2 <R2>

            What it is: Used for paired-end reads (e.g., from Illumina). -1 specifies the file with the forward reads, and -2 specifies the file with the reverse reads. This is the most common input type.

            Example: -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz

        -s <single>

            What it is: Used for single-end (unpaired) reads. Can be used alone or in addition to paired-end data.

        --trusted-contigs <fa>

            What it is: Provides a set of high-quality contigs in FASTA format to be used as a trusted "backbone" for the assembly.

            Critical Info: SPAdes will use these contigs to guide the assembly graph construction, aiming to extend them and bridge gaps using the read data. This is useful if you have a high-quality assembly of a very closely related reference genome.

        --untrusted-contigs <fa>

            What it is: Similar to the above, but SPAdes treats these contigs with more suspicion.

            Critical Info: If the read data strongly contradicts the structure of an untrusted contig, SPAdes has the authority to "break" it. This is useful for scaffolding with a more distant reference or a draft-quality assembly.

    Hybrid Assembly Options
    These flags are used for hybrid assembly, which combines the accuracy of short reads with the length of long reads to produce highly contiguous, often complete, genomes.

        --nanopore <fastq>

            What it is: Specifies an input file of long reads from an Oxford Nanopore sequencer.

        --pacbio <fastq>

            What it is: Specifies an input file of long reads from a Pacific Biosciences sequencer (this flag works for both older CLR and newer high-fidelity HiFi/CCS reads).

            Critical Info: Hybrid assembly is a powerful feature of SPAdes. It typically uses the long reads for scaffolding the de Bruijn graph built from the accurate short reads, resolving repeats and creating a much more complete final assembly.

    Pipeline Stage Control
    These flags allow you to run only specific parts of the SPAdes pipeline.

        --only-assembler

            What it is: Skips the initial read error-correction step (BayesHammer) and proceeds directly to the assembly stage.

            Critical Info: Use this if the reads are of very high quality or if you have already performed error correction using a different tool. This can save a significant amount of time.

        --only-error-correction

            What it is: Runs only the BayesHammer error-correction module and then exits.

            Critical Info: This is a very useful feature for pre-processing the reads. The corrected reads will be saved in the output directory and can then be used as input for SPAdes (with --only-assembler) or any other downstream analysis or different assembly tool.

For further details, please go to: https://ablab.github.io/spades/input.html
"""


# =====================================================================================
# main file structure is here
# =====================================================================================

def get_input_file():
    '''
    This function reads the command-line arguments and normalize the path.

    Why I used it:
      I wanted the script to be runnable without code edits. If the user did not provide a path,
      I defaulted to the current directory so they could just run the script in the right place.
      
    What I used:
      argparse.ArgumentParser to declare arguments and parse them,
      Path(...).expanduser().resolve() to normalize to an absolute path that worked well everywhere.
    '''
    parser = argparse.ArgumentParser(description="SPAdes input helper (simple structure)")
    parser.add_argument("path", nargs="?", help="Path to a file or a folder (folder recommended)")

    args = parser.parse_args()

    chosen_path = Path(args.path).expanduser().resolve() if args.path else Path(".").resolve()
    return chosen_path


def validate_input_file(file_path):
    '''
    I used this function to fail early with clear error messages.

    Why I used it:
      I wanted to separate “this path did not exist” from “this file was not readable,”
      because users fixed those problems differently.

    What I used:
      Path.exists() and Path.is_file() from pathlib for structure,
      os.access(path, os.R_OK) to check basic readability on the current system.
    '''
    p = Path(file_path)
    if not p.exists():
        print(f"[invalid] Path did not exist: {p}")
        return False
    if p.is_file() and not os.access(p, os.R_OK):
        print(f"[invalid] File was not readable: {p}")
        return False
    return True


def do_something_with_file(file_path):
    '''
    I used this function for the main flow:
      1) I cataloged files with previews and a quick classification.
      2) I asked the user to pick items by number and choose use / convert / skip.
      3) For FASTQ “use” files, I inferred pairs and estimated R1/R2 overlap.
      4) I printed a SPAdes template + short options list.
      5) I wrote a CSV report code for the catalog and pairing results.

    I wrote it this way to teach what each step achieved, while keeping the surface area small.
    '''
    # I resolved a folder to catalog. If the user pointed at a file, I cataloged its parent folder
    # so the user could still see neighboring files.
    folder = file_path if Path(file_path).is_dir() else Path(file_path).parent
    all_files = [f for f in sorted(folder.iterdir()) if f.is_file()]

    

    if not all_files:
        print("No files were found in this location.")
        return

    print(make_bold("[files found]"))
    catalog_rows = []
    for idx, fpath in enumerate(all_files, 1):
        ext = all_suffixes_as_one(fpath)
        kind, is_valid, why = classify_by_look(fpath)
        # I computed what the extension suggested so I could warn about mislabeled files.
        expected_kind = (
            "fastq" if ext in fastq_files else
            "fasta" if ext in fasta_files else
            {".bam": "bam", ".sam": "sam", ".cram": "cram",
             ".fast5": "hdf5", ".h5": "hdf5", ".bax.h5": "hdf5", ".sra": "sra"}.get(ext, "unknown")
        )
        extension_mismatch = (expected_kind != kind and kind != "unknown")

        print(make_bold(f"{idx:2d}. {fpath.name}"))
        print(f"    Extension : {ext}")
        print(f"    Detected  : {kind.upper()} ({'VALID' if is_valid else 'NOT VALID'})")
        if extension_mismatch:
            print(f"    {make_bold('[mismatch]')} extension suggested '{expected_kind}', but content looked like '{kind}'")
        print("    Preview   :")
        for line in show_preview(fpath).splitlines():
            print("      " + line)
        print("-" * 72)

        catalog_rows.append({
            "i": idx, "path": fpath, "kind": kind, "valid": is_valid,
            "why": why, "mismatch": extension_mismatch
        })

    # I asked people to select by number, because copying long paths was annoying and error-prone.
    raw_choice = input(make_bold("Select file numbers (e.g., 1,3,5-7 or all): ")).strip().lower()

    def parse_numbers(spec: str, max_n: int):
        '''
        I wrote this to turn input like "1,3,5-7" into [1,3,5,6,7].
        I accepted "all", "*" or "a" for convenience. I bounded indices to avoid crashes.
        '''
        if spec in ("all", "*", "a"):
            return list(range(1, max_n + 1))
        picked = set()
        for token in spec.split(","):
            token = token.strip()
            if not token:
                continue
            if "-" in token:
                a, b = token.split("-", 1)
                try:
                    a = int(a); b = int(b)
                    for k in range(a, b + 1):
                        if 1 <= k <= max_n:
                            picked.add(k)
                except:
                    pass
            else:
                try:
                    k = int(token)
                    if 1 <= k <= max_n:
                        picked.add(k)
                except:
                    pass
        return sorted(picked)

    picked_indices = parse_numbers(raw_choice, len(all_files))
    if not picked_indices:
        print("Nothing was selected.")
        return
    chosen_paths = [all_files[i - 1] for i in picked_indices]

    # I gave people three choices per file:
    #   use     → treat as ready for SPAdes (FASTQ as reads; FASTA as contigs)
    #   convert → print exact hints to get FASTQ first (BAM/CRAM/SAM/HDF5/SRA)
    #   skip    → ignore this file for now
    fastq_chosen, trusted_contigs, untrusted_contigs = [], [], []
    for one_path in chosen_paths:
        kind, is_valid, _ = classify_by_look(one_path)
        default_action = "use" if is_valid and kind in ("fastq", "fasta") else "convert"
        action = (input(f"{one_path.name} → Action [use/convert/skip] (default {default_action}): ")
                  .strip().lower() or default_action)
        if action == "skip":
            continue
        if action == "convert":
            # I printed conversion guidance. I did not actually run external tools in this script.
            if kind in ("bam", "sam", "cram"):
                extra = " (CRAM needed -T <ref.fa>)" if kind == "cram" else ""
                print(f"  To FASTQ{extra}:\n    samtools fastq {one_path} -1 reads_1.fastq -2 reads_2.fastq -s singles.fastq")
            elif kind == "sra":
                print(f"  To FASTQ:\n    fasterq-dump {one_path} -O output_folder")
            elif kind == "hdf5":
                base = one_path if one_path.is_dir() else one_path.parent
                print(f"  To FASTQ (Nanopore):\n    dorado basecaller dna_r9.4.1 -r {base}/ > nanopore.fastq")
            else:
                print("  (no built-in conversion hint)")
            continue
        if action == "use":
            if kind == "fastq":
                fastq_chosen.append(one_path)
            elif kind == "fasta":
                which = (input("Treat FASTA as [trusted/untrusted] (default untrusted): ")
                         .strip().lower() or "untrusted")
                (trusted_contigs if which == "trusted" else untrusted_contigs).append(one_path)

    # I grouped FASTQs by a “sample stub” and then tried to pair likely R1 with likely R2.
    # I also estimated overlap by comparing small sets of normalized read IDs.
    groups_by_sample = {}
    for fq in fastq_chosen:
        groups_by_sample.setdefault(infer_sample_stub(fq), []).append(fq)

    pairs, singles, pairing_rows = [], [], []
    used = set()
    for sample_id, group in groups_by_sample.items():
        r1_list = [p for p in group if decide_orientation(p) == "R1"]
        r2_list = [p for p in group if decide_orientation(p) == "R2"]
        while r1_list and r2_list:
            r1 = r1_list.pop(0); r2 = r2_list.pop(0)
            pairs.append((r1, r2)); used.update((r1, r2))
            a_ids = first_read_ids(r1, 100); b_ids = first_read_ids(r2, 100)
            inter = len(a_ids & b_ids)
            base = max(1, min(len(a_ids), len(b_ids)))
            overlap_pct = 100.0 * inter / base
            pairing_rows.append({
                "r1": r1.name, "r2": r2.name, "sample_id": sample_id,
                "overlap_pct": f"{overlap_pct:.1f}",
                "note": "ok" if overlap_pct >= 70 else "low overlap—check files"
            })
        for leftover in group:
            if leftover not in used:
                singles.append(leftover)

    # I printed only a template and an options cheat-sheet here on purpose, so users had to choose
    # sensible flags for their data rather than copy-pasting a fixed command.
    print("\nTEMPLATE:", make_bold(spades_template_line()))
    print("\n" + spades_options_block())

    # I wrote to the CSV so the user had a persistent summary they could email, archive, or compare between runs.
    report_csv = folder / "spades_input_report.csv"
    try:
        with open(report_csv, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(["index", "name", "path", "extension", "detected_kind", "valid", "mismatch", "note"])
            for row in catalog_rows:
                writer.writerow([
                    row["i"], row["path"].name, str(row["path"]),
                    all_suffixes_as_one(row["path"]),
                    row["kind"], "VALID" if row["valid"] else "NOT VALID",
                    "YES" if row["mismatch"] else "no", row["why"]
                ])
            writer.writerow([])
            writer.writerow(["r1", "r2", "sample_id", "overlap_pct", "note"])
            for pr in pairing_rows:
                writer.writerow([pr["r1"], pr["r2"], pr["sample_id"], pr["overlap_pct"], pr["note"]])
        print(f"Report was written to: {report_csv}")
    except Exception as err:
        print("I could not write the CSV report:", err)


# I used the standard Python main function call here to make the script runnable as
# a program and importable as a module without side effects.
if __name__ == "__main__":
    path_chosen = get_input_file()
    if validate_input_file(path_chosen):
        do_something_with_file(path_chosen)
