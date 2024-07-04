
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
Usage: KleTy.py [OPTIONS]

Options:
  -q, --query TEXT        query genome in fasta or fastq format. May be
                          gzipped.
  --ql TEXT               a list of query files. One query per line.
  -o, --prefix TEXT       prefix for output. Only work when there is only one
                          query. default: query filename
  -n, --n_proc INTEGER    number of process to use. default: 8
  -f, --plasmid_fragment  flag to predict plasmid fragment sharing < 50% with
                          the reference
  -m, --skip_gene         flag to skip AMR/VF searching. default: False
  -g, --skip_cgmlst       flag to skip cgMLST. default: False
  -p, --skip_plasmid      flag to skip plasmid typing. default: False
  --help                  Show this message and exit.
~~~~~~~~~~~~~~~~~

## Parameters:

| Parameter   | Explanation |
| ----------- | ----------- |
| **-q, --query** | Query genome. This can be in Fasta or Fastq format, and can be in plain text or GZIPped. |
| **--ql** | A list of query files. One query genome (file location) per line. KleTy will run these queries one by one and concatenate the outputs together. |
| **-o, --prefix** | Prefix for the outputs. There will be two files <prefix>.KleTy and <prefix>.cgMLST.profile.gz. Will use the prefix of the query file (or the ql file) if not specified.  |
| **-n, --n_proc** | Number of processes to use. Default: 8 |
| **-f, --plasmid_fragment** | Flag to predict less reliable plasmid fragments that share <50% (but >=30%) of the reference plasmid. |
| **-m, --skip_gene** | Flag to skip AMR/VF Searching. This step normaly taks ~ 15 seconds. |
| **-g, --skip_cgmlst** | Flag to skip cgMLST calling. This step normaly taks ~ 20 seconds.  |
| **-p, --skip_plasmid** | Flag to skip plasmid prediction. This step normaly taks ~ 30 seconds.  |


# Outputs:
## KleTy generates:

~~~~~~~~~~~~~
<prefix>.KleTy
~~~~~~~~~~~~~

### <prefix>.KleTy contains the genotyping results
~~~~~~~~~~~~~
$ cat CP015990.KleTy
INPUT   REPLICON        SPECIES HC1360.500.200.100.50.20.10.5.2 REFERENCE       PLASTYPE        COVERAGE        AMR:AMINOGLYCOSIDE      AMR:BETA-LACTAM AMR:CARBAPENEM  AMR:ESBL        AMR:INHIBITOR-RESISTANT AMR:COLISTIN    AMR:FOSFOMYCIN  AMR:MACROLIDE   AMR:PHENICOL    AMR:QUINOLONE   AMR:RIFAMYCIN   AMR:GLYCOPEPTIDES       AMR:SULFONAMIDE AMR:TETRACYCLINE        AMR:TIGECYCLINE AMR:TRIMETHOPRIM        AMR:BLA_INTRINSIC       STRESS:COPPER   STRESS:MERCURY  STRESS:NICKEL   STRESS:SILVER  STRESS:TELLURIUM STRESS:ARSENIC  STRESS:FLUORIDE STRESS:QUATERNARY_AMMONIUM      VIRULENCE:clb   VIRULENCE:iro   VIRULENCE:iuc   VIRULENCE:rmp   VIRULENCE:ybt   Others  REPLICON:INC_TYPE       REPLICON:MOB_TYPE       REPLICON:MPF_TYPE       ANNOTATION      CONTIGS
examples/CP015990.fna   ALL     Klebsiella_pneumoniae   10.10.10.10.ND.ND.ND.ND.ND      KLE_DA0156AA_AS -       -       aac(3)-IId^,aac(6')-Ib-cr.v2^,aadA16*   OXA-1   KPC-2   -       -       -       -       mphA    catB3.v2        GyrA-83F,GyrA-87A,ParC-80I,qnrA3^       arr-3   -       sul1    -       -       dfrA27  SHV-28^ -       merA,merE,merR_Ps,merT  -       -       -       -       -       qacEdelta1      -       -       -       -       fyuA_26,irp1_275,irp2_30,ybtA_78,ybtE_58,ybtP_75,ybtQ_88,ybtS_115,ybtT_26,ybtU_129,ybtX_73      -       IncR    -       MPF_T   -       -
examples/CP015990.fna   P1      -       -       CP059309.1      PT_361,PC_361   84.9    aac(6')-Ib-cr.v2^,aadA16*       OXA-1   KPC-2   -       -       -       -       mphA    catB3.v2        qnrA3^  arr-3   -       sul1    -       -       dfrA27 --       merA,merE,merR_Ps,merT  -       -       -       -       -       qacEdelta1      -       -       -       -       -       -       IncR    -       -       Klebsiella_pneumoniae_strain_Kp46596_plasmid_pKp46596-3,_complete_sequence      CP015991.1
examples/CP015990.fna   Others  -       -       -       -       -       aac(3)-IId^     -       KPC-2   -       -       -       -       mphA    -       GyrA-83F,GyrA-87A,ParC-80I      -       -       -       -       -       -       SHV-28^ -      --       -       -       -       -       -       -       -       -       -       fyuA_26,irp1_275,irp2_30,ybtA_78,ybtE_58,ybtP_75,ybtQ_88,ybtS_115,ybtT_26,ybtU_129,ybtX_73      -       -       -       MPF_T   -       -
~~~~~~~~~~~~~

The columns are:
* **INPUT**: Filename of the input. Used to recognize query assemblies. 
* **REPLICON**: Type of the replicon. It can be:
~~~~~~~~~~~~~
  ALL: A summary of the query.
  P<n>: One plasmid per row. Will not be reported with '-p'. 
  Others: Summary of the AMR/VF genes that are not in plasmids (likely carried by the chromosome). Will not be reported with '-p'. 
~~~~~~~~~~~~~
  
* **SPECIES**: Species designation of the query, inferred based on its cgMLST profile. Will not be reported with '-g'. 
* **HC1360.500.200.100.50.20.10.5.2**: HierCC cluster designation of the query based on the cgMLST profile. HC1360 approximately equals to clonal complex (CC) in MLST. Lower HC levels were used for sub-population clusterings. Numbers after HC indicate the criteria of the single-linkage clustering. Will not be reported with '-g'. 
* **REFERENCE**: Accession of the reference for predicted plasmid. Will not be reported with '-p'. 
* **PLASTYPE**: PT (plasmid type) and PC (plasmid cluster) of the predicted plasmid. Will not be reported with '-p'. 
* **COVERAGE**: Coverage of the plasmid to the reference. Will not be reported with '-p'. 
* **AMR:AMINOGLYCOSIDE**: Predicted genes/mutations encoding resistance to AMINOGLYCOSIDE. 
* **AMR:BETA-LACTAM**: Predicted genes/mutations encoding resistance to BETA-LACTAM. 
* **AMR:CARBAPENEM**: Predicted genes/mutations encoding resistance to CARBAPENEM. 
* **AMR:ESBL**: Predicted genes/mutations encoding Extended-spectrum beta-lactamases (ESBLs). 
* **AMR:INHIBITOR-RESISTANT**: Predicted genes/mutations encoding resistance to Beta-Lactamase inhibitors. 
* **AMR:COLISTIN**: Predicted genes/mutations encoding resistance to COLISTIN. 
* **AMR:FOSFOMYCIN**: Predicted genes/mutations encoding resistance to FOSFOMYCIN. 
* **AMR:MACROLIDE**: Predicted genes/mutations encoding resistance to MACROLIDE. 
* **AMR:PHENICOL**: Predicted genes/mutations encoding resistance to PHENICOL. 
* **AMR:QUINOLONE**: Predicted genes/mutations encoding resistance to QUINOLONE. 
* **AMR:RIFAMYCIN**: Predicted genes/mutations encoding resistance to RIFAMYCIN. 
* **AMR:GLYCOPEPTIDES**: Predicted genes/mutations encoding resistance to GLYCOPEPTIDES. 
* **AMR:SULFONAMIDE**: Predicted genes/mutations encoding resistance to SULFONAMIDE. 
* **AMR:TETRACYCLINE**: Predicted genes/mutations encoding resistance to TETRACYCLINE. 
* **AMR:TIGECYCLINE**: Predicted genes/mutations encoding resistance to TIGECYCLINE. 
* **AMR:TRIMETHOPRIM**: Predicted genes/mutations encoding resistance to TRIMETHOPRIM. 
* **AMR:BLA_INTRINSIC**: Predicted intrinsic beta-lactamase in Klebsiella. 
* **STRESS:COPPER**: Predicted genes encoding resistance to COPPER. 
* **STRESS:MERCURY**: Predicted genes encoding resistance to MERCURY. 
* **STRESS:NICKEL**: Predicted genes encoding resistance to NICKEL. 
* **STRESS:SILVER**: Predicted genes encoding resistance to SILVER. 
* **STRESS:TELLURIUM**: Predicted genes encoding resistance to TELLURIUM. 
* **STRESS:ARSENIC**: Predicted genes encoding resistance to ARSENIC. 
* **STRESS:FLUORIDE**: Predicted genes encoding resistance to FLUORIDE. 
* **STRESS:QUATERNARY_AMMONIUM**: Predicted genes encoding resistance to QUATERNARY_AMMONIUM. 
* **VIRULENCE:clb**: colibactin (clb)
* **VIRULENCE:iro**: salmochelin (iro)
* **VIRULENCE:iuc**: aerobactin (iuc)
* **VIRULENCE:rmp**: hypermucoidy (rmpA, rmpA2)
* **VIRULENCE:ybt**: yersiniabactin (ybt)
* **Others**: Other resistances. 
* **REPLICON:INC_TYPE**: INC type of the plasmid. 
* **REPLICON:MOB_TYPE**: MOB type of the plasmid. 
* **REPLICON:MPF_TYPE**: MPF type of the plasmid. 
* **ANNOTATION**: Annotations of the predicted plasmids. 
* **CONTIGS**: Contigs associated with the predicted plasmids. 


| Column      | Explanation |
| ----------- | ----------- |
| **INPUT**   | Filename of the input. Used to recognize query assemblies |
| **REPLICON** | Type of the replicon. It can be: 
~~~~~~~~~~~~~
  ALL: A summary of the query.
  P<n>: One plasmid per row. Will not be reported with '-p'. 
  Others: Summary of the AMR/VF genes that are not in plasmids (likely carried by the chromosome). Will not be reported with '-p'. 
~~~~~~~~~~~~~
|
| **SPECIES** | Species designation of the query, inferred based on its cgMLST profile. Will not be reported with '-g'. |
| **HC1360.500.200.100.50.20.10.5.2** | HierCC cluster designation of the query based on the cgMLST profile. HC1360 approximately equals to clonal complex (CC) in MLST. Lower HC levels were used for sub-population clusterings. Numbers after HC indicate the criteria of the single-linkage clustering. Will not be reported with '-g'. |
| **REFERENCE** | Accession of the reference for predicted plasmid. Will not be reported with '-p'. |
| **PLASTYPE** | PT (plasmid type) and PC (plasmid cluster) of the predicted plasmid. Will not be reported with '-p'. |
| **COVERAGE** | Coverage of the plasmid to the reference. Will not be reported with '-p'. |
| **AMR:AMINOGLYCOSIDE** | Predicted genes/mutations encoding resistance to AMINOGLYCOSIDE. |
| **AMR:BETA-LACTAM** | Predicted genes/mutations encoding resistance to BETA-LACTAM. |
| **AMR:CARBAPENEM** | Predicted genes/mutations encoding resistance to CARBAPENEM. |
| **AMR:ESBL** | Predicted genes/mutations encoding Extended-spectrum beta-lactamases (ESBLs). |
| **AMR:INHIBITOR-RESISTANT** | Predicted genes/mutations encoding resistance to Beta-Lactamase inhibitors. |
| **AMR:COLISTIN** | Predicted genes/mutations encoding resistance to COLISTIN. |
| **AMR:FOSFOMYCIN** | Predicted genes/mutations encoding resistance to FOSFOMYCIN. |
| **AMR:MACROLIDE** | Predicted genes/mutations encoding resistance to MACROLIDE. |
| **AMR:PHENICOL** | Predicted genes/mutations encoding resistance to PHENICOL. |
| **AMR:QUINOLONE** | Predicted genes/mutations encoding resistance to QUINOLONE. |
| **AMR:RIFAMYCIN** | Predicted genes/mutations encoding resistance to RIFAMYCIN. |
| **AMR:GLYCOPEPTIDES** | Predicted genes/mutations encoding resistance to GLYCOPEPTIDES. |
| **AMR:SULFONAMIDE** | Predicted genes/mutations encoding resistance to SULFONAMIDE. |
| **AMR:TETRACYCLINE** | Predicted genes/mutations encoding resistance to TETRACYCLINE. |
| **AMR:TIGECYCLINE** | Predicted genes/mutations encoding resistance to TIGECYCLINE. |
| **AMR:TRIMETHOPRIM** | Predicted genes/mutations encoding resistance to TRIMETHOPRIM. |
| **AMR:BLA_INTRINSIC** | Predicted intrinsic beta-lactamase in Klebsiella. |
| **STRESS:COPPER** | Predicted genes encoding resistance to COPPER. |
| **STRESS:MERCURY** | Predicted genes encoding resistance to MERCURY. |
| **STRESS:NICKEL** | Predicted genes encoding resistance to NICKEL. |
| **STRESS:SILVER** | Predicted genes encoding resistance to SILVER. |
| **STRESS:TELLURIUM** | Predicted genes encoding resistance to TELLURIUM. |
| **STRESS:ARSENIC** | Predicted genes encoding resistance to ARSENIC. |
| **STRESS:FLUORIDE** | Predicted genes encoding resistance to FLUORIDE. |
| **STRESS:QUATERNARY_AMMONIUM** | Predicted genes encoding resistance to QUATERNARY_AMMONIUM. |
| **VIRULENCE:clb** | colibactin (clb) |
| **VIRULENCE:iro** | salmochelin (iro) |
| **VIRULENCE:iuc** | aerobactin (iuc) |
| **VIRULENCE:rmp** | hypermucoidy (rmpA, rmpA2) |
| **VIRULENCE:ybt** | yersiniabactin (ybt) |
| **Others** | Other resistances |
| **REPLICON:INC_TYPE** | INC type of the plasmid. |
| **REPLICON:MOB_TYPE** | MOB type of the plasmid. |
| **REPLICON:MPF_TYPE** | MPF type of the plasmid. |
| **ANNOTATION** | Annotations of the predicted plasmids. |
| **CONTIGS** | Contigs associated with the predicted plasmids. |


# Citation and Reproduction Instructions

### Reproduction Instructions
All data required for reproduction of the analysis were distributed in this repository under
https://github.com/zheminzhou/KleTy/tree/main/db


These includes:
* plasmids.repr.clu.gz - **IMPORTANT**. A mapping table that specifies correlations between plasmids and PT/PCs. 
* HierCC.tsv.gz - A tab-delimited table consisting of HierCC results for all ~70,000 genomes
* klebsiella.cgmlst - A list of core genes used in the dcgMLST scheme
* klebsiella.refsets.fas.gz - reference alleles for all pan genes (for calling new alleles)
* klebsiella.species - A mapping table that specifies correlations between genomes and Klebsiella species
* profile.parq - Allelic profiles of all ~70,000 genomes in parquet format, and can be read using the Pandas library (https://pandas.pydata.org/docs/reference/api/pandas.read_parquet.html). 
* stress_CDS.gz - reference sequences for resistance to metal/biocides
* traditional_lasmid_type.fas.gz - reference sequences for INC/MOB/MPF types of the plasmids. 
* kleborate/* - reference sequences from kleborate. 


