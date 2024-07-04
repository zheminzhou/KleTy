
# KleTy (Klebsiella typer for the core genome and plasmids)
KleTy is a tool to type Klebsiella genome assemblies for: 
* core genome MLST (cgMLST) for detailed genotyping of the core genome
* Hierarchical clusters (HierCC) that represents natural population
* Plasmid prediction and classification (PC)
* hypervirulence associated loci
* antimicrobial resistance determinants


# INSTALLATION:

KleTy was developed and tested in Python >=3.8. It depends on several Python libraries: 
~~~~~~~~~~
click
numba
numpy
pandas
biopython
pyarrow
fastparquet
~~~~~~~~~~

All libraries can be installed using pip: 

~~~~~~~~~~
pip install click numba numpy pandas biopython pyarrow fastparquet
~~~~~~~~~~
KleTy also calls NCBI-BLAST+:

~~~~~~~~~~
ncbi-blast+
~~~~~~~~~~

Which can be installed via 'apt' in UBUNTU:
~~~~~~~~~~
sudo apt install -y ncbi-blast+
~~~~~~~~~~

The whole environment can also be installed in conda:


~~~~~~~~~~
conda create --name dty python==3.11
conda activate dty
conda install -c conda-forge biopython numba numpy pandas click pyarrow fastparquet
conda install -c bio-conda blast
~~~~~~~~~~

The installation process normally finishes in <10 minutes. 

NOTE: Please make sure that "makeblastdb" and "blastn" are all in the PATH environment variable (can be run without pointing to their actual location). 

When run for the first time, KleTy will automatically download the reference plasmids from https://zenodo.org/records/12590507/files/plasmids.repr.fas.gz
This will only run once. But note that the file is fairly large (816 MB), and will take a long time to download. 

Alternatively, for those who have difficulty downloading the file within the pipeline. Please download the file by yourself, and copy it into "db/" under the KleTy folder. Then run

~~~~~~~~~~
gzip -d plasmids.repr.fas.gz
makeblastdb -in plasmids.repr.fas -dbtype nucl
~~~~~~~~~~

to generate the required database. 



# Quick Start (with examples)
## Get allelic and HierCC callings
~~~~~~~~~~~
$ cd /path/to/KleTy/
$ python KleTy.py -q examples/CP015990.fna
~~~~~~~~~~~

The whole calculation finishes in ~1 minutes with 8 CPU threads (~2.5 minutes with one CPU thread). 
The screen output will be like:

~~~~~~~~~~~
07/04/2024 04:44:56 AM Running query: examples/CP015990.fna
07/04/2024 04:44:56 AM  Searching VF/STRESS genes...
07/04/2024 04:45:06 AM  Done.
07/04/2024 04:45:06 AM  Searching AMR genes...
07/04/2024 04:45:16 AM  Done.
07/04/2024 04:45:16 AM  Searching plasmids...
07/04/2024 04:45:51 AM  Done.
07/04/2024 04:45:51 AM  Running cgMLST...
07/04/2024 04:46:09 AM  Done.
~~~~~~~~~~~

And there are two outputs (see below for explanation):
~~~~~~~~~~~
CP015990.KleTy
CP015990.cgMLST.profile.gz
~~~~~~~~~~~


# USAGE:
## KleTy.py - allelic callings and HierCC clusters & species predictions

~~~~~~~~~~~~~~
$ python KleTy.py --help
Usage: KleTy.py [OPTIONS]

Options:
  -q, --query TEXT      query genome in fasta or fastq format. May be gzipped.
  --ql TEXT             a list of query files. One query per line.
  -o, --prefix TEXT     prefix for output. Only work when there is only one
                        query. default: query filename
  -n, --n_proc INTEGER  number of process to use. default: 8
  -m, --skip_gene       flag to skip AMR/VF searching. default: False
  -g, --skip_cgmlst     flag to skip cgMLST. default: False
  -p, --skip_plasmid    flag to skip plasmid typing. default: False
  --help                Show this message and exit.
~~~~~~~~~~~~~~~~~

# Outputs:
## KleTy generates:

~~~~~~~~~~~~~
<prefix_of_input>.KleTy
~~~~~~~~~~~~~

### <prefix_of_input>.KleTy contains the genotyping results
~~~~~~~~~~~~~
$ cat CP015990.KleTy
REPLICON        SPECIES HC1360.500.200.100.50.20.10.5.2 REFERENCE       PLASTYPE        COVERAGE        AMR:AMINOGLYCOSIDE      AMR:AMIKACIN    AMR:APRAMYCIN   AMR:GENTAMICIN  AMR:HYGROMYCIN  AMR:KANAMYCIN   AMR:SPECTINOMYCIN       AMR:STREPTOMYCIN        AMR:TOBRAMYCIN AMR:BETA-LACTAM  AMR:CARBAPENEM  AMR:CEPHALOSPORIN       AMR:ESBL        AMR:INHIBITOR-RESISTANT AMR:COLISTIN    AMR:FOSFOMYCIN  AMR:LINCOSAMIDE AMR:MACROLIDE   AMR:PHENICOL    AMR:CHLORAMPHENICOL     AMR:FLORFENICOL AMR:KIRROMYCIN  AMR:PULVOMYCIN  AMR:QUINOLONE   AMR:RIFAMYCIN   AMR:STREPTOTHRICIN      AMR:SULFONAMIDE AMR:TETRACYCLINE        AMR:TIGECYCLINE AMR:TRIMETHOPRIM        AMR:BLEOMYCIN   STRESS:COPPER   STRESS:MERCURY  STRESS:NICKEL   STRESS:SILVER   STRESS:TELLURIUM        STRESS:ARSENIC  STRESS:FLUORIDE STRESS:QUATERNARY_AMMONIUM      VIRULENCE:clb   VIRULENCE:iro   VIRULENCE:iuc   VIRULENCE:rmp   VIRULENCE:ybt   Others  REPLICON:INC_TYPE       REPLICON:MOB_TYPE       REPLICON:MPF_TYPE       ANNOTATION      CONTIGS
CP015990:ALL    Klebsiella_pneumoniae   10.10.10.10.ND.ND.ND.ND.ND      KLE_DA0156AA_AS -       -       aac(3)-IId,aac(6')-Ib-cr5,aadA16        aac(6')-Ib-cr5  -       aac(3)-IId      -       aac(6')-Ib-cr5  -       aadA16  aac(6')-Ib-cr5  blaKPC-2,blaOXA-1,blaSHV-28    blaKPC-2 blaOXA-1        -       -       -       -       -       mph(A)  catB3   catB3   -       -       -       aac(6')-Ib-cr5,qnrA3    arr-3,rpoB_V146F        -       sul1    -       -       dfrA27  -       silA,silR       MerP_Gneg,merA,merD,merE,merR_Ps,merR_Ps(*Premature),merT       -       silA,silR       -       arsB_pKW301     -       qacEdelta1      -       -       -       -       fyuA_26,irp1_275,irp2_30,ybtA_73,ybtE_58,ybtP_75,ybtQ_88,ybtS_115,ybtT_26,ybtU_129,ybtX_73      acrF,emrD,kdeA(AMR:EFFLUX)      IncR    -       MPF_T  --
CP015990:P1     -       -       CP059309.1      PT_361,PC_361   84.9    aac(6')-Ib-cr5,aadA16   aac(6')-Ib-cr5  -       -       -       aac(6')-Ib-cr5  -       aadA16  aac(6')-Ib-cr5  blaKPC-2,blaOXA-1       blaKPC-2        blaOXA-1        -       -       -       -      -mph(A)  catB3   catB3   -       -       -       aac(6')-Ib-cr5,qnrA3    arr-3   -       sul1    -       -       dfrA27  -       -       MerP_Gneg,merA,merD,merE,merR_Ps,merR_Ps(*Premature),merT       -       -       -       -       -       qacEdelta1      -       -      --       -       -       IncR    -       -       Klebsiella_pneumoniae_strain_Kp46596_plasmid_pKp46596-3,_complete_sequence      CP015991.1
CP015990:Others -       -       -       -       -       aac(3)-IId      -       -       aac(3)-IId      -       -       -       -       -       blaKPC-2,blaSHV-28      blaKPC-2        -       -       -       -       -       -       mph(A)  -       -       -       -      --       rpoB_V146F      -       -       -       -       -       -       silA,silR       -       -       silA,silR       -       arsB_pKW301     -       -       -       -       -       -       fyuA_26,irp1_275,irp2_30,ybtA_73,ybtE_58,ybtP_75,ybtQ_88,ybtS_115,ybtT_26,ybtU_129,ybtX_73      acrF,emrD,kdeA(AMR:EFFLUX)      -       -       MPF_T   -       -
~~~~~~~~~~~~~


# Citation and Reproduction Instructions

### Reproduction Instructions
All data required for reproduction of the analysis were distributed in this repository under
https://github.com/zheminzhou/KleTy/tree/main/db


These includes:
* klebsiella.refsets.fas.gz - reference alleles for all pan genes (for calling new alleles)
* klebsiella.cgmlst - A list of core genes used in the dcgMLST scheme
* profile.parq - Allelic profiles of all ~70,000 genomes in parquet format, and can be read using the Pandas library (https://pandas.pydata.org/docs/reference/api/pandas.read_parquet.html). 
* HierCC.tsv.gz - A tab-delimited table consisting of HierCC results for all ~70,000 genomes
* klebsiella.species - A mapping table that specifies correlations between genomes and Klebsiella species
* virulence.fasta - reference sequences for hypervirulence genes as used in kleborate
* AMR* - reference sequences for AMR genes and mutations as used in AMRfinderPlus
* resfams_db - A mapping file that gives more accurate prediction of ESBLs
* plasmids/* - Reference sequences for INC and MOB types. 


