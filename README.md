# CasCollect

**CasCollect** is a pipeline for the detection of Cas genes and CRISPR arrays from unassembled raw illumina sequencing reads.
This pipeline relies on python and perl as well as requiring the installation of several pieces of software within the users path:

+ BBTools
+ Seqtk
+ FragGeneScan
+ HMMER
+ VSEARCH
+ SPAdes
+ CRISPRCasFinder

## optional arguments:

  -h, --help          show this help message and exit

### flags for read trimming/cleaning
```
  --trim              runs quality trim, remove adapters, merge overlapping
                      reads
  --clean             runs quality trim, remove adapters and reference
                      matching reads [use: -ref file.fasta for user defined
                      reference file], merge overlapping reads
  -ref file.fasta     user defined contaminant reference genome/sequences
                      [default: human genome]
```
### arguments for read files
```
  -fwd file.fastq     fastq file of forward reads [use with -rev]
  -rev file.fastq     fastq file of reverse reads [use with -fwd]
  -single file.fastq  fastq file of unpaired reads [use in place of -fwd and
                      -rev]
```
### flags and arguments for seed file(s) generation
```
  --noprot            disable search for Cas or other protein genes
                      withinunassembled reads for assembly [use with -hmm for
                      user-defined hmm file]
  -hmm hmm dir        directory containing user defined protein hmm file(s)
                      for protein search within unassembled reads [default:
                      Cas protein hmm]
  --nucl              search for CRISPR or other nucleotide sequences
                      withinunassembled reads for assembly [must use with
                      -query for user-defined fasta file]
  -query file.fasta   user defined fasta DNA file of sequences
  --seed              user-defined set of seeds for read subset expansion
  -define file.fasta  user-defined fasta DNA file of seed
```
### arguments for read subset expansion
```
  -cycle number       number of cycles for seed expansion with vsearch
  -match number       percent match for vsearch
```
### flags for read assembly
```
  --noassembly        disable assembly of reads extracted from protein and/or
                      nucleotide searches [disables downstream annotation]
  --meta              run assembly for metagenomic data
```
### flags for annotation
```
  --noannotate        disable CRISPR arrays and Cas protein gene annotation of
                      assembled reads from protein and/or nucleotide searches
```
### arguments for computational resources
```
  -cpu threads        #CPUs for all tasks
  -mem RAM            Gb of RAM for all tasks
```
### arguments for output location
```
  -out folder         path to and folder name for all outputs
```
### flags for program version
```
  --version, -v       show program's version number and exit
```

## Minimal arguments for running CasCollect
### With paired-end data is:
```
python CasCollect.py -fwd file.fastq -rev file.fastq -out folder
```
### With single-end data is:
```
python CasCollect.py - single file.fastq -out folder
```
This will run **CasCollect** for the specified fastq files to search the default 120 Cas protein HMMs with 5 cycles of seed expansion on 1 cpu and 20 Gb RAM outputting all resulting into the specified folder.
