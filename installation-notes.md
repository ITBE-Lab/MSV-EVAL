# Installation Notes
The below commands install all required software for reproducing the experiments in "State-of-the-art structural variant calling: What went conceptually wrong and how to fix it?". As environment, we suggest an Ubuntu 20.04.1 installation.


## Basics

    sudo apt-get -y install build-essential git cmake python3 python3-dev python3-pip zlib1g zlib1g-dev autoconf clang libc++-dev libc++abi-dev


## Install PostgreSQL

    sudo sh -c 'echo "deb http://apt.postgresql.org/pub/repos/apt $(lsb_release -cs)-pgdg main" > /etc/apt/sources.list.d/pgdg.list'
    wget --quiet -O - https://www.postgresql.org/media/keys/ACCC4CF8.asc | sudo apt-key add -
    sudo apt-get update
    sudo apt-get -y install postgresql-12 libpq-dev postgresql-server-dev-12


## Install PostGIS

    sudo apt install postgis postgresql-12-postgis-3
    sudo -u postgres psql
    CREATE EXTENSION postgis;
    \q


## Configure Postgres User

    # log in without password
    sudo -u postgres psql 
    # set password to admin
    \password postgres
in the password prompt enter “admin” as password (without the double quotes)

    # quit
    \q

    # allow md5 connections
    sudo nano /etc/postgresql/12/main/pg_hba.conf

using the editor, replace the line "local all postgres peer" with "local all postgres md5"

    # restart service
    sudo service postgresql restart

## (Optional) Install pgAdmin4

    pip install pgadmin4

## Install MSV as a Python Module

    export PYTHONPATH=$PYTHONPATH:~/buildMA
if you intend to install MSV permanently, you should add the above export statement to your login bash-script (~/.bashrc file). If you omit this addition to your login bash-script, the above export directive gets lost after shell closure.

## Dowload & Compile MSV

    git clone https://github.com/ITBE-Lab/MA.git
    cd MA
    git checkout svCaller
    cd ..
    mkdir buildMA
    cd buildMA
    CC="clang" CXX="clang++" cmake -DWITH_PYTHON=ON ../MA/
    make -j 8
    cd ..


## Install MSV-EVAL

    git clone https://github.com/ITBE-Lab/MSV-EVAL.git
    sudo pip3 install bokeh==1.4.0


## Install Minimap2

    git clone https://github.com/lh3/minimap2.git
    cd minimap2
    make -j 8
    cd ..

## Install NGMLR

    wget https://github.com/philres/ngmlr/releases/download/v0.2.7/ngmlr-0.2.7-linux-x86_64.tar.gz
    tar xvzf ngmlr-0.2.7-linux-x86_64.tar.gz
    mv ngmlr-0.2.7 ngmlr


## Install SNIFFLES

    git clone https://github.com/fritzsedlazeck/Sniffles.git
    cd Sniffles
    mkdir build
    cd build
    cmake ..
    make -j 8
    mv ../bin/sniffles-core-* ../bin/sniffles-core
    cd ../..


## Install Delly

    mkdir delly
    cd delly
    wget https://github.com/dellytools/delly/releases/download/v0.8.6/delly_v0.8.6_linux_x86_64bit -O delly
    chmod +x delly
    cd ..

## Install GraphAligner and VG

    conda install -c bioconda graphaligner
    wget https://github.com/vgteam/vg/releases/download/v1.36.0/vg
    chmod +x vg


## Install Gridss

    conda create -n gridss gridss
    conda activate gridss


## Install Manta

    conda install manta


## Install Bcftools

    sudo apt-get -y install bcftools


## Install Samtools

    sudo apt-get -y install samtools


## Install SURVIVOR

    git clone https://github.com/fritzsedlazeck/SURVIVOR.git
    cd SURVIVOR/Debug
    make -j 8
    cd ../..


## Install DWGSIM

    git clone --recursive https://github.com/nh13/DWGSIM.git
    cd DWGSIM
    make -j 8
    cd ..


## Setup Folder Structure

    sudo mkdir /MAdata
    sudo chown msv:msv /MAdata
    mkdir -p /MAdata/genome/yeasts/UFRJ50816/fasta
    mkdir -p /MAdata/genome/yeasts/YPS138/fasta
    mkdir -p /MAdata/genome/human/GRCh38.p12/fasta
    mkdir -p /MAdata/genome/reconstructed/yeast/UFRJ50816/ma
    mkdir -p /MAdata/sv_caller_analysis/svs_hidden_to_aligners/reads
    mkdir -p /MAdata/sv_caller_analysis/svs_hidden_to_aligners/sam
    mkdir -p /MAdata/sv_caller_analysis/ambiguities_of_atomic_sv/sam
    mkdir -p /MAdata/sv_caller_analysis/ambiguities_of_atomic_sv/vcf
    mkdir -p /MAdata/sv_caller_analysis/yeast_analysis
    mkdir -p /MAdata/sv_caller_analysis/yeast_analysis/gridss
    mkdir -p /MAdata/ena/simulated/UFRJ50816/Illumina-250
    mkdir -p /MAdata/ena/simulated/UFRJ50816/pacbio_CCS
    mkdir -p /MAdata/ena/simulated/UFRJ50816/oxfNano
    mkdir -p /MAdata/tmp
alternatively, you can reconfigure the script MS-EVAL/sv_util/settings.py 

## Download Yeast Genomes & Build Indices
The below statements download the yeast genomes required for reproducing Fig. 4 and Table 1 of the manuscript.

    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/079/115/GCA_002079115.1_ASM207911v1/GCA_002079115.1_ASM207911v1_genomic.fna.gz -O /MAdata/genome/yeasts/YPS138/fasta/genome.fna.gz
    gunzip /MAdata/genome/yeasts/YPS138/fasta/genome.fna.gz

    mkdir /MAdata/genome/yeasts/YPS138/minimap
    ~/minimap2/minimap2 -x map-pb -d /MAdata/genome/yeasts/YPS138/minimap/genome.map-pb.mmi /MAdata/genome/yeasts/YPS138/fasta/genome.fna
    mkdir /MAdata/genome/yeasts/YPS138/ngmlr
    cp /MAdata/genome/yeasts/YPS138/fasta/genome.fna /MAdata/genome/yeasts/YPS138/ngmlr/genome.fna
    mkdir /MAdata/genome/yeasts/YPS138/ma
    ~/buildMA/maCMD --Create_Index /MAdata/genome/yeasts/YPS138/fasta/genome.fna,/MAdata/genome/yeasts/YPS138/ma,genome

    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/079/145/GCA_002079145.1_ASM207914v1/GCA_002079145.1_ASM207914v1_genomic.fna.gz -O /MAdata/genome/yeasts/UFRJ50816/fasta/genome.fna.gz
    gunzip /MAdata/genome/yeasts/UFRJ50816/fasta/genome.fna.gz

    mkdir /MAdata/genome/yeasts/UFRJ50816/minimap
    ~/minimap2/minimap2 -x map-pb -d /MAdata/genome/yeasts/UFRJ50816/minimap/genome.map-pb.mmi /MAdata/genome/yeasts/UFRJ50816/fasta/genome.fna
    mkdir /MAdata/genome/yeasts/UFRJ50816/ngmlr
    cp /MAdata/genome/yeasts/UFRJ50816/fasta/genome.fna /MAdata/genome/yeasts/UFRJ50816/ngmlr/genome.fna
    mkdir /MAdata/genome/yeasts/UFRJ50816/ma
    ~/buildMA/maCMD --Create_Index /MAdata/genome/yeasts/UFRJ50816/fasta/genome.fna,/MAdata/genome/yeasts/UFRJ50816/ma,genome


## Download Human Genome & Build Indices
The below statements download the human genomee required for reproducing Fig. 1 and Fig. 2 of the manuscript.
The downloading and index building can take several hours.

    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.27_GRCh38.p12/GCA_000001405.27_GRCh38.p12_genomic.fna.gz -O /MAdata/genome/human/GRCh38.p12/fasta/genome.fna.gz
    gunzip /MAdata/genome/human/GRCh38.p12/fasta/genome.fna.gz

    mkdir /MAdata/genome/human/GRCh38.p12/minimap
    ~/minimap2/minimap2 -x map-pb -d /MAdata/genome/human/GRCh38.p12/minimap/genome.map-pb.mmi /MAdata/genome/human/GRCh38.p12/fasta/genome.fna
    ~/minimap2/minimap2 -x asm10 -d /MAdata/genome/human/GRCh38.p12/minimap/genome.asm10.mmi /MAdata/genome/human/GRCh38.p12/fasta/genome.fna
    mkdir /MAdata/genome/human/GRCh38.p12/ngmlr
    cp /MAdata/genome/human/GRCh38.p12/fasta/genome.fna /MAdata/genome/human/GRCh38.p12/ngmlr/genome.fna
    mkdir /MAdata/genome/human/GRCh38.p12/ma
    ~/buildMA/maCMD --Create_Index /MAdata/genome/human/GRCh38.p12/fasta/genome.fna,/MAdata/genome/human/GRCh38.p12/ma,genome

    mkdir /MAdata/genome/human/GRCh38.p12/vg
    printf "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" > /MAdata/genome/human/GRCh38.p12/vg/empty.vcf
    vg construct -r /MAdata/genome/human/GRCh38.p12/fasta/genome.fna -v /MAdata/genome/human/GRCh38.p12/vg/empty.vcf > /MAdata/genome/human/GRCh38.p12/vg/genome.vg

## adjust RAM max usage of KSW

The data for Table 1 of the manuscript are computed via Dynamic Programming. By default this computation is disabled, since it writes large files (>250GB) to the disk. For enabling it, the following modifications of the script ~/MSV-Eval/sv_util/settings.py are requried:
- Set the variable run_ksw to True.
- Set the variable ksw_file_system_min_gb_size to at least 4GB less than the RAM of your machine (larger values increase comparison speed). E.g. if your machine has 32GB RAM, the maximal (and recommended) value is 28GB.


## Run Experiments

    cd ~/MSV-EVAL

    # Fig. 1 (requires previous download of the human genome)
    python3 ambiguities_of_atomic_sv/main.py

    # Fig. 2 (requires previous download of the human genome)
    python3 svs_hidden_to_aligners/main.py

    # Fig. 4 & Table 1 (requires previous download of the yeast genomes)
    python3 yeast/main.py


