#!/usr/bin/env python3
import argparse
import os
import subprocess
import time
import glob
from distutils.spawn import find_executable

# start timing for files
start_files = time.time()

parser = argparse.ArgumentParser(description='Cas protein read collection, targeted assembly, and annotation '
                                             '(CasCollect) is a flexible pipeline designed for the detection of Cas '
                                             'protein genes and CRISPR direct-repeat from unassembled next-generation '
                                             'short-read sequencing data.'
                                             'Requires:'
                                             'BBTools'
                                             'Seqtk'
                                             'FragGeneScan'
                                             'HMMER'
                                             'VSEARCH'
                                             'SPAdes'
                                             'CRISPRCasFinder')
# flags for read trimming/cleaning
parser.add_argument('--trim', action='store_true', help='runs quality trim, remove adapters, merge overlapping reads')
parser.add_argument('--clean', action='store_true', help='runs quality trim, remove adapters and reference matching'
                                                         ' reads [use: -ref file.fasta for user defined reference'
                                                         ' file], merge overlapping reads')
parser.add_argument('-ref', metavar='file.fasta', type=str, help='user defined contaminant reference genome/sequences'
                                                                 ' [default: human genome]')
# arguments for read files
parser.add_argument('-fwd', metavar='file.fastq', type=str, help='fastq file of forward reads [use with -rev]')
parser.add_argument('-rev', metavar='file.fastq', type=str, help='fastq file of reverse reads [use with -fwd]')
parser.add_argument('-single', metavar='file.fastq', type=str, help='fastq file of unpaired reads [use in place of -fwd'
                                                                    ' and -rev]')
# arguments for seed file(s) generation
parser.add_argument('--noprot', action='store_true', help='disable search for Cas or other protein genes within'
                                                          'unassembled reads for assembly [use with -hmm for user'
                                                          '-defined hmm file]')
parser.add_argument('-hmm', metavar='hmm dir', type=str, help='directory containing user defined protein hmm file(s)'
                                                              ' for protein search within unassembled reads [default:'
                                                              ' Cas protein hmm]')
parser.add_argument('--nucl', action='store_true', help='search for CRISPR or other nucleotide sequences within'
                                                        'unassembled reads for assembly [must use with -query for user'
                                                        '-defined fasta file]')
parser.add_argument('-query', metavar='file.fasta', type=str, help='user defined fasta DNA file of sequences')
parser.add_argument('--seed', action='store_true', help='user-defined set of seeds for read subset expansion')
parser.add_argument('-define', metavar='file.fasta', type=str, help='user-defined fasta DNA file of seed')
# read subset expansion
parser.add_argument('-cycle', metavar='number', type=int, default=5, help='number of cycles for seed expansion '
                                                                          'with vsearch')
parser.add_argument('-match', metavar='number', type=str, default=str(0.95), help='percent match for vsearch')
# arguments for read assembly
parser.add_argument('--noassembly', action='store_true', help='disable assembly of reads extracted from protein and/or'
                                                              ' nucleotide searches [disables downstream annotation]')
parser.add_argument('--meta', action='store_true', help='run assembly for metagenomic data')
# arguments for annotation
parser.add_argument('--noannotate', action='store_true', help='disable CRISPR arrays and Cas protein gene annotation of'
                                                              ' assembled reads from protein and/or nucleotide'
                                                              ' searches')
# arguments for computational resources
parser.add_argument('-cpu', metavar='threads', type=str, default=str(1), help='#CPUs for all tasks')
parser.add_argument('-mem', metavar='RAM', type=int, default=20, help='Gb of RAM for all tasks')
# arguments for output location
parser.add_argument('-out', metavar='folder', type=str, required=True, help='path to and folder name for all outputs')
# flags for program version
parser.add_argument('--version', '-v', action='version', version='%(prog)s 1.0')

args = parser.parse_args()

trim = args.trim
clean = args.clean
ref = args.ref

fwd = args.fwd
rev = args.rev
single = args.single

noprot = args.noprot
hmm = args.hmm
nucl = args.nucl
query = args.query
seed = args.seed
define = args.define

cycle = args.cycle
match = args.match

noassembly = args.noassembly
meta = args.meta
noannotate = args.noannotate

cpu = args.cpu
mem = str(args.mem*1000)

out = args.out

# absolute path to files and scripts
program_path = (os.path.dirname(os.path.realpath(__file__)))
adapters = os.path.dirname(find_executable('bbduk.sh')) + '/resources/adapters.fa'
ccfsum = program_path + '/bin/ccfSum.pl'

# check for paired and unpaired fastq files
if str(fwd) == 'None' and str(rev) == 'None' and str(single) == 'None':
    # no fastq files provided, kill program
    print('no paired files or unpaired fastq file provided')
    exit()
elif str(fwd) != 'None' and str(single) != 'None':
    # paired and unpaired reads as inout
    print('choose only paired or unpaired reads as input')
    exit()
elif str(fwd) != 'None' and str(rev) == 'None':
    # no reverse fastq files provided for paired, kill program
    print('reverse fastq file needed')
    exit()
elif str(fwd) == 'None' and str(rev) != 'None':
    # no forward fastq files provided for paired, kill program
    print('forward fastq file needed')
    exit()
elif str(fwd) != 'None' and str(rev) != 'None':
    print('paired-end reads as input')
    read_input = 'paired'
    fwd_file = os.path.abspath(fwd)
    rev_file = os.path.abspath(rev)
elif str(single) != 'None':
    print('unpaired reads as input')
    read_input = 'single'
    single_file = os.path.abspath(single)
else:
    print('issues with input files')
    exit()

# ref for cleaning
if str(ref) != 'None':
    ref = os.path.abspath(ref)

# hmm file for Cas protein genes or other proteins
if str(hmm) == 'None':
    hmm = program_path + '/ref'
else:
    hmm = os.path.abspath(hmm)

# query file CRISPR or other nucleotide sequences
if nucl and str(query) != 'None':
    query = os.path.abspath(query)
elif nucl and str(query) == 'None':
    print('failed to detect query file')
    exit()
else:
    print('no nucleotide query')

# user-defined seed file
if seed and str(define) != 'None':
    query = os.path.abspath(define)
elif seed and str(define) == 'None':
    print('failed to detect user-defined seed file')
    exit()
else:
    print('no user-defined seed file')

# check if folder exists, append folder name
if os.path.exists(out):
    print('folder already exists in this location, appending folder name')
    out = out + '_' + str(time.time())
else:
    out = out
os.mkdir(out)
os.chdir(out)

# trimming or cleaning command
if read_input == 'paired' and (trim or clean):
    os.mkdir('fastq')
    if trim:
        print('trimming paired-end reads')
        ref_input = 'ref=' + adapters
    elif clean:
        print('cleaning paired-end reads')
        ref_input = 'ref=' + adapters + ',' + ref
    subprocess.call(['bbduk.sh', 'in1=' + fwd_file,  'in2=' + rev_file, 'out1=fastq/fwd_file.fastq',
                     'out2=fastq/rev_file.fastq', 'outs=fastq/single.fastq', ref_input, 'ktrim=r',
                     'minlen=25', 'minlenfraction=0.333', 'mink=11', 'overwrite=true', 'k=23', 'threads=' + cpu])
    subprocess.call(['bbmerge.sh', 'in1=fastq/fwd_file.fastq', 'in2=fastq/rev_file.fastq',
                     'out=fastq/merged.fastq', 'outu1=fastq/fwd_file_unmerged.fastq',
                     'outu2=fastq/rev_file_unmerged.fastq', 'threads=' + cpu])
    single_fastq = open('fastq/single_file.fastq', 'w')
    subprocess.call(['cat', 'fastq/single.fastq', 'fastq/merged.fastq'], stdout=single_fastq)
elif read_input == 'single' and (trim or clean):
    os.mkdir('fastq')
    if trim:
        print('trimming paired-end reads')
        ref_input = 'ref=' + adapters
    elif clean:
        print('cleaning paired-end reads')
        ref_input = 'ref=' + adapters + ',' + ref
    subprocess.call(['bbduk.sh', 'in=' + single_file, 'out=fastq/single_file.fastq', ref_input,
                     'ktrim=r', 'minlen=25', 'minlenfraction=0.333', 'mink=11', 'overwrite=true', 'k=23', 'threads=' + cpu])
else:
    print('no trimming or cleaning requested')

# convert fastq to fasta and combine, need for protein or nucleotide seed generation
os.mkdir('fasta')
reads = open('fasta/reads.fasta', 'w')
# for paired-end reads trimmed or cleaned
if read_input == 'paired' and (trim or clean):
    forward = open('fasta/fwd_file_unmerged.fasta', 'w')
    subprocess.call(['seqtk', 'seq', '-a', 'fastq/fwd_file_unmerged.fastq'], stdout=forward)
    reverse = open('fasta/rev_file_unmerged.fasta', 'w')
    subprocess.call(['seqtk', 'seq', '-a', 'fastq/rev_file_unmerged.fastq'], stdout=reverse)
    singleton = open('fasta/single_file.fasta', 'w')
    subprocess.call(['seqtk', 'seq', '-a', 'fastq/single_file.fastq'], stdout=singleton)
    subprocess.call(['cat', 'fasta/fwd_file_unmerged.fasta', 'fasta/rev_file_unmerged.fasta',
                     'fasta/single_file.fasta'], stdout=reads)
# for unpaired reads trimmed or cleaned
elif read_input == 'single' and (trim or clean):
    subprocess.call(['seqtk', 'seq', '-a', 'fastq/single_file.fastq'], stdout=reads)
# for paired-end reads not trimmed or cleaned
elif read_input == 'paired' and not trim and not clean:
    forward = open('fasta/fwd_file.fasta', 'w')
    subprocess.call(['seqtk', 'seq', '-a', fwd_file], stdout=forward)
    reverse = open('fasta/rev_file.fasta', 'w')
    subprocess.call(['seqtk', 'seq', '-a', rev_file], stdout=reverse)
    subprocess.call(['cat', 'fasta/fwd_file.fasta', 'fasta/rev_file.fasta'], stdout=reads)
# for unpaired reads not trimmed or cleaned
elif read_input == 'single' and not trim and not clean:
    subprocess.call(['seqtk', 'seq', '-a', single_file], stdout=reads)
else:
    print('failed to convert fastq to fasta')
    exit()

# end timing for files
end_files = time.time()

# start timing for seed generation
start_seed_gen = time.time()

# protein seed generation command
if not noprot:
    # translate fasta
    train = os.path.dirname(find_executable('FGS+'))
    os.mkdir('prot_seeds')
    subprocess.call(['FGS+', '-s', 'fasta/reads.fasta', '-o', 'prot_seeds/translated', '-w', '0', '-r', train +
                     '/train/', '-t', 'illumina_1', '-p', cpu, '-m', mem])
    # search with hmm
    FNULL = open(os.devnull, 'w')
    hmm_list = glob.glob(hmm + '/*.hmm')
    for each_hmm in hmm_list:
        tblout = each_hmm.replace(hmm, 'prot_seeds/')
        tblout = tblout.replace('.hmm', '')
        subprocess.call(['hmmsearch', '--tblout', tblout, '--noali', '--cpu', cpu, each_hmm,
                         'prot_seeds/translated.faa'], stdout=FNULL)
        prot_tblout_1 = open(tblout + '.1', 'w')
        subprocess.call(['cut', '-d', '_', '-f', '1', tblout], stdout=prot_tblout_1)
        prot_tblout_2 = open(tblout, 'w')
        subprocess.call(['tail', '-n', '+4', tblout + '.1'], stdout=prot_tblout_2)
        prot_tblout_3 = open(tblout + '.1', 'w')
        subprocess.call(['head', '-n', '-10', tblout], stdout=prot_tblout_3)
        os.remove(tblout)
    tblout = glob.glob('prot_seeds/*.1')
    with open('prot_seeds/prot_hits', 'w') as outfile:
        for reads in tblout:
            with open(reads, "r") as infile:
                outfile.write(infile.read())
    uniq_prot_hits = open('prot_seeds/uniq_prot_hits', 'w')
    subprocess.call(['sort', '-u', 'prot_seeds/prot_hits'], stdout=uniq_prot_hits)
    # extract fasta reads from protein seed generation
    prot_hits = open('prot_seeds/prot_hits.fasta', 'w')
    # for paired-end trimmed or cleaned reads
    if read_input == 'paired' and (trim or clean):
        prot_hits_fwd = open('prot_seeds/prot_hits_fwd.fasta', 'w')
        subprocess.call(['seqtk', 'subseq', 'fasta/fwd_file_unmerged.fasta', 'prot_seeds/uniq_prot_hits'],
                        stdout=prot_hits_fwd)
        prot_hits_rev = open('prot_seeds/prot_hits_rev.fasta', 'w')
        subprocess.call(['seqtk', 'subseq', 'fasta/rev_file_unmerged.fasta', 'prot_seeds/uniq_prot_hits'],
                        stdout=prot_hits_rev)
        prot_hits_single = open('prot_seeds/prot_hits_single.fasta', 'w')
        subprocess.call(['seqtk', 'subseq', 'fasta/single_file.fasta', 'prot_seeds/uniq_prot_hits'],
                        stdout=prot_hits_single)
        subprocess.call(['cat', 'prot_seeds/prot_hits_fwd.fasta', 'prot_seeds/prot_hits_rev.fasta',
                        'prot_seeds/prot_hits_merged.fasta', 'prot_seeds/prot_hits_single.fasta'], stdout=prot_hits)
    # for paired-end reads not trimmed or cleaned
    elif read_input == 'paired' and not trim and not clean:
        prot_hits_fwd = open('prot_seeds/prot_hits_fwd.fasta', 'w')
        subprocess.call(['seqtk', 'subseq', 'fasta/fwd_file.fasta', 'prot_seeds/uniq_prot_hits'], stdout=prot_hits_fwd)
        prot_hits_rev = open('prot_seeds/prot_hits_rev.fasta', 'w')
        subprocess.call(['seqtk', 'subseq', 'fasta/rev_file.fasta', 'prot_seeds/uniq_prot_hits'], stdout=prot_hits_rev)
        subprocess.call(['cat', 'prot_seeds/prot_hits_fwd.fasta', 'prot_seeds/prot_hits_rev.fasta'], stdout=prot_hits)
    # for unpaired reads
    elif read_input == 'single':
        subprocess.call(['seqtk', 'subseq', 'fasta/reads.fasta', 'prot_seeds/uniq_prot_hits'], stdout=prot_hits)
    else:
        print('failed protein seed generation')
        exit()
else:
    print('no protein seed generation requested')

# nucleotide seed generation command
if nucl:
    os.mkdir('nucl_seeds')
    subprocess.call(['vsearch', '--usearch_global', query, '--db', 'fasta/reads.fasta', '--minseqlength', '1', '--id',
                     '1', '--strand', 'both', '--maxaccepts', '100', '--uc_allhits', '--userfields', 'target',
                     '--threads', cpu, '--userout', 'nucl_seeds/hits'])
    hits_uniq = open('nucl_seeds/hits.uniq', 'w')
    subprocess.call(['sort', '-u', 'nucl_seeds/hits'], stdout=hits_uniq)
    nucl_hits = open('nucl_seeds/nucl_hits.fasta', 'w')
    subprocess.call(['seqtk', 'subseq', 'fasta/reads.fasta', 'nucl_seeds/hits.uniq'], stdout=nucl_hits)
else:
    print('no nucleotide seed generation requested')

# assign protein and nucleotide seeds
# protein and nucleotide and user-defined seeds
if not noprot and nucl and seed:
    combined_seeds = open('fasta/seeds.fasta', 'w')
    subprocess.call(['cat', 'prot_seeds/prot_hits.fasta', 'nucl_seeds/nucl_hits.fasta', define], stdout=combined_seeds)
    seeds = 'fasta/seeds.fasta'
# protein and user-defined seeds
elif not noprot and not nucl and seed:
    combined_seeds = open('fasta/seeds.fasta', 'w')
    subprocess.call(['cat', 'prot_seeds/prot_hits.fasta', define], stdout=combined_seeds)
    seeds = 'fasta/seeds.fasta'
# nucleotide and user-defined seeds
elif noprot and nucl and seed:
    combined_seeds = open('fasta/seeds.fasta', 'w')
    subprocess.call(['cat', 'nucl_seeds/nucl_hits.fasta', define], stdout=combined_seeds)
    seeds = 'fasta/seeds.fasta'
# protein and nucleotide seeds
elif not noprot and nucl and not seed:
    combined_seeds = open('fasta/seeds.fasta', 'w')
    subprocess.call(['cat', 'prot_seeds/prot_hits.fasta', 'nucl_seeds/nucl_hits.fasta'], stdout=combined_seeds)
    seeds = 'fasta/seeds.fasta'
# protein seeds
elif not noprot and not nucl and not seed:
    seeds = 'prot_seeds/prot_hits.fasta'
# nucleotide seeds
elif noprot and nucl and not seed:
    seeds = 'nucl_seeds/nucl_hits.fasta'
# user-defined seeds
elif noprot and not nucl and seed:
    seeds = define
else:
    print('seed determination failed')
    exit()


# cycles of seed expansion
db = 'fasta/reads.fasta'

for y in range(1,cycle + 1):
    subprocess.call(['vsearch', '--usearch_global', seeds, '--db', db, '--minseqlength', '1', '--id',
                     match, '--strand', 'both', '--maxaccepts', '1000000', '--uc_allhits', '--dbmatched',
                     'fasta/seeds_' + str(y) + '.fasta', '--dbnotmatched', 'fasta/remainder_' + str(y) + '.fasta',
                     '--userfields', 'target', '--userout', 'fasta/seeds_' + str(y), '--threads', cpu])
    seeds_uniq = open('fasta/seeds_' + str(y) + '.uniq', 'w')
    subprocess.call(['sort', '-u', 'fasta/seeds_' + str(y)], stdout=seeds_uniq)
    seeds = 'fasta/seeds_' + str(y) + '.fasta'
    db = 'fasta/remainder_' + str(y) + '.fasta'

seed_headers = glob.glob('fasta/seeds_*.uniq')
if not noprot:
    seed_headers.append('prot_seeds/uniq_prot_hits')
with open('fasta/seed_headers', 'w') as outfile:
    for x in seed_headers:
        with open(x) as infile:
            outfile.write(infile.read())
seed_headers_uniq = open('fasta/seed_headers.uniq', 'w')
subprocess.call(['sort', '-u', 'fasta/seed_headers'], stdout=seed_headers_uniq)

# extract fastq read seeds
if read_input == 'paired':
    if trim or clean:
        fwd_fastq = 'fastq/fwd_file_unmerged.fastq'
        rev_fastq = 'fastq/rev_file_unmerged.fastq'
        single_fastq = 'fastq/single_file.fastq'
    elif not trim and not clean:
        fwd_fastq = fwd
        rev_fastq = rev
    fwd_fastq_file = open('fastq/fwd_seeds.fastq', 'w')
    subprocess.call(['seqtk', 'subseq', fwd_fastq, 'fasta/seed_headers.uniq'], stdout=fwd_fastq_file)
    rev_fastq_file = open('fastq/rev_seeds.fastq', 'w')
    subprocess.call(['seqtk', 'subseq', rev_fastq, 'fasta/seed_headers.uniq'], stdout=rev_fastq_file)
elif read_input == 'single':
    if trim or clean:
        single_fastq = 'fastq/single_file.fastq'
    elif not trim and not clean:
        single_fastq = single
if read_input == 'single' or (read_input == 'paired' and (trim or clean)):
    single_fastq_file = open('fastq/single_seeds.fastq', 'w')
    subprocess.call(['seqtk', 'subseq', single_fastq, 'fasta/seed_headers.uniq'], stdout=single_fastq_file)

# end timing for seed generation
end_seed_gen = time.time()

# start timing for assembly
start_assembly = time.time()

# assembly
if not noassembly:
    os.mkdir('assembly')
    # for paired-end trimmed or cleaned reads
    if read_input == 'paired':
        if trim or clean:
            fwd_reads = 'fastq/fwd_seeds.fastq'
            rev_reads = 'fastq/rev_seeds.fastq'
            single_reads = 'fastq/single_seeds.fastq'
            if meta:
                subprocess.call(['spades.py', '--meta', '-1', fwd_reads, '-2', rev_reads, '-s', single_reads, '-t',
                                 cpu, '-m', str(args.mem), '-o', 'assembly'])
            elif not meta:
                subprocess.call(['spades.py', '-1', fwd_reads, '-2', rev_reads, '-s', single_reads, '-t', cpu, '-m',
                                 str(args.mem), '-o', 'assembly'])
            else:
                print('issue with assembly command')
        elif not trim and not clean:
            fwd_reads = 'fastq/fwd_seeds.fastq'
            rev_reads = 'fastq/rev_seeds.fastq'
            if meta:
                subprocess.call(['spades.py', '--meta', '-1', fwd_reads, '-2', rev_reads, '-t', cpu, '-m',
                                 str(args.mem), '-o', 'assembly'])
            elif not meta:
                subprocess.call(['spades.py', '-1', fwd_reads, '-2', rev_reads, '-t', cpu, '-m', str(args.mem), '-o',
                                 'assembly'])
            else:
                print('issue with assembly command')
        else:
            print('failed to determine input files for assembly')
    elif read_input == 'single':
        single_reads = 'fastq/single_seeds.fastq'
        if meta:
            subprocess.call(['spades.py', '--meta', '-s', single_reads, '-t', cpu, '-m', str(args.mem), '-o',
                             'assembly'])
        elif not meta:
            subprocess.call(['spades.py', '-s', single_reads, '-t', cpu, '-m', str(args.mem), '-o', 'assembly'])
        else:
            print('issue with assembly command')
    else:
        print('assembly failed')
        exit()
else:
    print('no assembly requested')

# end timing for assembly
end_assembly = time.time()

# start timing for annotation
start_annotation = time.time()

# CRISPR and cas protein gene annotation
if noassembly or noannotate:
    print('no annotation requested')
elif not noassembly and not noannotate:
    subprocess.call(['perl', ccfsum, 'assembly/scaffolds.fasta', 'search', 'spades', cpu])
else:
    print('annotation failed')
    exit()

# end timing for annotation
end_annotation = time.time()


print('end of code')

# write timings
os.chdir('..')
timing_file = open('timing.' + out, 'w')
timing_file.write('trim=' + str(trim) + '\n'
                  + 'clean=' + str(clean) + '\n'
                  + 'ref=' + str(ref) + '\n'
                  + 'fwd=' + str(fwd) + '\n'
                  + 'rev=' + str(rev) + '\n'
                  + 'single=' + str(single) + '\n'
                  + 'noprot=' + str(noprot) + '\n'
                  + 'hmm=' + str(hmm) + '\n'
                  + 'nucl=' + str(nucl) + '\n'
                  + 'query=' + str(query) + '\n'
                  + 'seed=' + str(seed) + '\n'
                  + 'define=' + str(define) + '\n'
                  + 'cycle=' + str(cycle) + '\n'
                  + 'noassembly=' + str(noassembly) + '\n'
                  + 'noannotate=' + str(noannotate) + '\n'
                  + 'cpu=' + cpu + '\n'
                  + 'mem=' + str(args.mem) + ' Gb\n'
                  + 'out=' + out + '\n'
                  + '\n' + '----' + '\n' + '\n'
                  + 'time for files' + '\n'
                  + str(end_files - start_files) + '\n'
                  + 'time for seed generation' + '\n'
                  + str(end_seed_gen - start_seed_gen) + '\n'
                  + 'time for assembly' + '\n'
                  + str(end_assembly - start_assembly) + '\n'
                  + 'time for annotation' + '\n'
                  + str(end_annotation - start_annotation) + '\n'
                  + '\n' + '----' + '\n' + '\n'
                  + 'time for CasCollect' + '\n'
                  + str(end_annotation - start_files))
