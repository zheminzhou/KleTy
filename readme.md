
# KleTy (Klebsiella typer for the core genome and plasmids)
Typing engine for core genome and plasmids in Klebsiella


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
KleTy also calls two 3rd party programs:

~~~~~~~~~~
ncbi-blast+
diamond
~~~~~~~~~~

Both can be installed via 'apt' in UBUNTU:
~~~~~~~~~~
sudo apt install -y ncbi-blast+ diamond
~~~~~~~~~~

The whole environment can also be installed in conda:


~~~~~~~~~~
conda create --name dty python==3.11
conda activate dty
conda install -c conda-forge biopython numba numpy pandas click pyarrow fastparquet
conda install -c bio-conda blast diamond
~~~~~~~~~~

The installation process normally finishes in <10 minutes. 

NOTE: Please make sure that "makeblastdb", "blastn", and "diamond" are all in the PATH environment variable (can be run without pointing to their actual location). 

Finally, format the plasmid reference database:
~~~~~~~~~~~
$ cd /path/to/KleTy/
$ cd db
$ gzip -d plasmids.repr.fas.gz
$ makeblastdb -in plasmids.repr.fas -dbtype nucl
~~~~~~~~~~~






# Quick Start (with examples)
## Get allelic and HierCC callings
~~~~~~~~~~~
$ cd /path/to/KleTy/
$ python KleTy.py -q examples/CP015990.fna
~~~~~~~~~~~

The whole calculation finishes in <3 minutes with 8 CPU threads. 


# USAGE:
## KleTy.py - allelic callings and HierCC clusters & species predictions

~~~~~~~~~~~~~~
$ Usage: KleTy.py [OPTIONS]

Options:
  -q, --query TEXT      query genome in fasta or fastq format. May be gzipped.
                        [required]
  -o, --prefix TEXT     prefix for output. default: query filename
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


