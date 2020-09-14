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

Before running, **CasCollect** will check for the necessary software in the users path and denote any missing software.  Alternatively, the **Check.py** script can be run to check for the necessary software in the users path with any missing software automatically downloaded and extracted.  The user is then responsible for installing each piece of software and ensuring it is in their path.

If useful, please cite: [Podlevsky JD, Hudson CM, Jerilyn A. Timlin JA, Williams KP. 2020. CasCollect: Targeted Assembly of CRISPR-Associated Operons from High-throughput Sequencing Data. NAR Genom Bioinform 2: lqaa063.](https://academic.oup.com/nargab/article/2/3/lqaa063/5901064)

## optional arguments:
```
  -h, --help          show this help message and exit
```
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

The default run for **CasCollect** is in protein mode.

## Protein Mode:

### With paired-end data is:
```
python CasCollect.py -fwd file.fastq -rev file.fastq -out folder
```
### With single-end data is:
```
python CasCollect.py -single file.fastq -out folder
```
This will run **CasCollect** for the specified fastq files to search the default 120 Cas protein HMMs with 5 cycles of seed expansion on 1 cpu and 20 Gb RAM outputting all resulting into the specified folder.

### To specify the HMM file(s):
```
python CasCollect.py -fwd file.fastq -rev file.fastq -out folder -hmm hmm_dir
```
This will run **CasCollect** using all the HMM files in the specified directory.

## DNA Mode:

### With paired-end data is:
```
python CasCollect.py -fwd file.fastq -rev file.fastq -out folder --noprot --nucl -query file.fasta
```
### With single-end data is:
```
python CasCollect.py -single file.fastq -out folder --noprot --nucl -query file.fasta
```
This will run **CasCollect** only for the DNA mode, searching with the query file.  No protein search.

## User-defined Mode:

### With paired-end data is:
```
python CasCollect.py -fwd file.fastq -rev file.fastq -out folder --noprot --seed -define file.fasta
```
### With single-end data is:
```
python CasCollect.py -single file.fastq -out folder --noprot --seed -define file.fasta
```
This will run **CasCollect** only for the User-defined mode, using the define file as the seeds for seed exapnsion and assembly.  No protein or DNA search.

## Combined Protein/DNA/User-defined mode:

### With paired-end data is:
```
python CasCollect.py -fwd file.fastq -rev file.fastq -out folder --nucl -query file.fasta --seed -define file.fasta
```
### With single-end data is:
```
python CasCollect.py -single file.fastq -out folder --nucl -query file.fasta --seed -define file.fasta
```
This will run **CasCollect** in Protein, DNA, and User-defined mode.  All the seeds generated from this mode will be combined and used for the downstream seed expansion and assembly.
