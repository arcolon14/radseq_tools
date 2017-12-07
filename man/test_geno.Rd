\name{test_geno}
\alias{test_geno}
\docType{data}
\title{
Test Genomic Sequence
}
\description{
gzipped, multi-line FASTA file containing the first 100Kb of six linkage groups of the Three Spine Stickleback (Gasterosteus aculeatus) reference genome.
}
\usage{./path/to/test_geno.fa.gz}
\format{
\code{>groupI dna:group group:BROADS1:groupI:1:100000:1}
\code{\code{CTCTTTGTTTTCAGGTGTGGAATGTGCTTTCTACCACGGCTACAAATACTACAAAGGATG}}
\code{TAGTAATGCTACTATATAGAACTAATACTGCAAATACTAAAACTACCACTAAAACTACTA}
\code{...}

\code{>groupII dna:group group:BROADS1:groupII:1:100000:1}
\code{TGACTATTAAGGCGTTTGCAGGCTGAGAGAAGCCAGTCTTGAATGCTACCCCCTTTTGAG}
\code{GAGAGATCCAGAATGACTGTTGAGACACTTGAATCAAATTCATCTTGGTCACTGAGCAGC}
\code{...}


\code{>groupIV dna:group group:BROADS1:groupIV:1:100000:1}
\code{GCTAGCCTGTTTAAACCAAACCATCGGTGTGTATGATTACTTGCGCCCACACCCGGTCTC}
\code{TGACGTCTGGGATTCTGGCAGCGGCGAGTGTATGTTCCCCTCCCCTCACGGGGTGCTGTA}
\code{...}

\code{...}


}

\source{
http://www.ensembl.org/Gasterosteus_aculeatus/Info/Index
}

\examples{
# Define path to test FASTA (test_geno.fa.gz)
path_fasta <- '/path/to/test_geno.fa.gz'

# Feed it into `process_fasta`
seq <- process_fasta(path_fasta)

