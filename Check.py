#!/usr/bin/env python3
import shutil
import wget
import subprocess
import tarfile
import zipfile
import os

print('CasCollect requires BBTools, Seqtk, FragGeneScanPlus, HMMER, VSEARCH, SPAdes, and CRISPRCasFinder.\nChecking if these programs are in your PATH.')

bbduk = shutil.which('bbduk.sh')
seqtk = shutil.which('seqtk')
FGS = shutil.which('FGS+')
hmmsearch = shutil.which('hmmsearch')
vsearch = shutil.which('vsearch')
spades = shutil.which('spades.py')
crisprcasfinder = shutil.which('CRISPRCasFinder.pl')

print('\nChecking for BBTools:')
if bbduk != None:
    print('BBTools is installed.')
else:
    print('BBTools is not installed, will attempt to download.\nPlease finish installation.')
    try:
        wget.download('https://downloads.sourceforge.net/project/bbmap/BBMap_38.84.tar.gz')
        bbduk = 1
    except:
        print('Cannot automatically download BBTools.\nPlease manually download from <https://downloads.sourceforge.net/project/bbmap/BBMap_38.84.tar.gz>')
        bbduk = 0
if bbduk == 1:
    bbtools = tarfile.open('BBMap_38.84.tar.gz')
    bbtools.extractall('.')
    bbtools.close()
    os.remove('BBMap_38.84.tar.gz')

print('\nChecking for Seqtk:')
if seqtk != None:
    print('Seqtk is installed.')
else:
    print('Seqtk is not installed, will attempt to download.\nPlease finish installation.')
    try:
        wget.download('https://github.com/lh3/seqtk/archive/master.zip')
        seqtk = 1
    except:
        print('Cannot automatically download Seqtk.\nPlease manually download from <https://github.com/lh3/seqtk/archive/master.zip>')
        seqtk = 0
if seqtk == 1:
    seqtk = zipfile.ZipFile('seqtk-master.zip')
    seqtk.extractall('.')
    seqtk.close()
    os.remove('seqtk-master.zip')

print('\nChecking for FragGeneScanPlus:')
if FGS != None:
    print('FragGeneScanPlus is installed.')
else:
    print('FragGeneScanPlus is not installed, will attempt to download.\nPlease finish installation.')
    try:
        wget.download('https://github.com/hallamlab/FragGeneScanPlus/archive/master.zip')
        FGS = 1
    except:
        print('Cannot automatically download FragGeneScanPlus.\nPlease manually download from <https://github.com/hallamlab/FragGeneScanPlus/archive/master.zip>')
        FGS = 0
if FGS == 1:
    FGS = zipfile.ZipFile('FragGeneScanPlus-master.zip')
    FGS.extractall('.')
    FGS.close()
    os.remove('FragGeneScanPlus-master.zip')

print('\nChecking for HMMER:')
if hmmsearch != None:
    print('HMMER is installed.')
else:
    print('HMMER is not installed, will attempt to download.\nPlease finish installation.')
    try:
        wget.download('http://eddylab.org/software/hmmer/hmmer.tar.gz')
        hmmsearch = 1
    except:
        print('Cannot automatically download HMMER.\nPlease manually download from <http://eddylab.org/software/hmmer/hmmer.tar.gz>')
        hmmsearch = 0
if hmmsearch == 1:
    hmmsearch = tarfile.open('hmmer.tar.gz')
    hmmsearch.extractall('.')
    hmmsearch.close()
    os.remove('hmmer.tar.gz')

print('\nChecking for VSEARCH:')
if vsearch != None:
    print('VSEARCH is installed.')
else:
    print('VSEARCH is not installed, will attempt to download.\nPlease finish installation.')
    try:
        wget.download('https://github.com/torognes/vsearch/archive/master.zip')
        vsearch = 1
    except:
        print('Cannot automatically download VSEARCH.\nPlease manually download from <https://github.com/torognes/vsearch/archive/master.zip>')
        vsearch = 0
if vsearch == 1:
    vsearch = zipfile.ZipFile('vsearch-master.zip')
    vsearch.extractall('.')
    vsearch.close()
    os.remove('vsearch-master.zip')

print('\nChecking for SPAdes:')
if spades != None:
    print('SPAdes is installed.')
else:
    print('SPAdes is not installed, will attempt to download.\nPlease finish installation.')
    try:
        wget.download('https://github.com/ablab/spades/archive/spades_3.14.1.zip')
        spades = 1
    except:
        print('Cannot automatically download SPAdes.\nPlease manually download from <https://github.com/ablab/spades/archive/spades_3.14.1.zip>')
        spades = 0
if spades == 1:
    spades = zipfile.ZipFile('spades-spades_3.14.1.zip')
    spades.extractall('.')
    spades.close()
    os.remove('spades-spades_3.14.1.zip')

print('\nChecking for CRISPRCasFinder:')
if crisprcasfinder != None:
    print('CRISPRCasFinder is installed.')
else:
    print('CRISPRCasFinder is not installed, will attempt to download.\nPlease finish installation.')
    try:
        wget.download('https://crisprcas.i2bc.paris-saclay.fr/Home/DownloadFile?filename=CRISPRCasFinder.zip')
        crisprcasfinder = 1
    except:
        print('Cannot automatically download CRISPRCasFinder.\nPlease manually download from <https://crisprcas.i2bc.paris-saclay.fr/Home/DownloadFile?filename=CRISPRCasFinder.zip>')
        crisprcasfinder = 0
if crisprcasfinder == 1:
    crisprcasfinder = zipfile.ZipFile('crisprcasfinder-DownloadFile?filename=CRISPRCasFinder.zip')
    crisprcasfinder.extractall('.')
    crisprcasfinder.close()
    os.remove('crisprcasfinder-DownloadFile?filename=CRISPRCasFinder.zip')

print('\nend of code')