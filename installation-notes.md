Virtual Box Version 6.1.16
Ubuntu 20.04.1
10 GB disk 16GB ram
User: MSV Password: default

## Basics

    sudo apt-get install build-essential git cmake python3 python3-dev python3-pip zlib1g zlib1g-dev autoconf clang libc++-dev libc++abi-dev


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
    admin
    admin
    # quit
    \q

    # allow md5 connections
    sudo nano /etc/postgresql/12/main/pg_hba.conf
    >
    > replace line 
    > local all postgres peer
    > with
    > local all postgres md5
    >

    # restart service
    sudo service postgresql restart


## Install MSV as Python Module

    # add this line to your ~/.bashrc file to permanently install MSV
    export PYTHONPATH=$PYTHONPATH:~/buildMA


## Dowload & Compile MSV

    git clone @todo
    git checkout svCaller
    mkdir buildMA
    cd buildMA
    CC="clang" CXX="clang++" cmake -DWITH_PYTHON=ON ../MA/
    make -j 32
    cd ..


## iInstall MSV-EVAL

    git clone @todo
    @todo add survivor sample file & adjust path...
    sudo pip3 install bokeh==1.4.0


## Install Minimap2

    git clone https://github.com/lh3/minimap2.git
    cd minimap2
    make -j 32
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
    make -j 32
    mv sniffles-core-* sniffles-core
    cd ../..


## Install Delly

    mkdir delly
    cd delly
    wget https://github.com/dellytools/delly/releases/download/v0.8.6/delly_v0.8.6_linux_x86_64bit -O delly
    chmod +x delly
    cd ..


## Install bcftools

    sudo apt-get install bcftools


## Install samtools

    sudo apt-get install samtools


## Install SURVIVOR

    git clone https://github.com/fritzsedlazeck/SURVIVOR.git
    cd SURVIVOR/Debug
    make -j 32
    cd ../..


## Install DWGSIM

    git clone --recursive https://github.com/nh13/DWGSIM.git
    cd DWGSIM
    make -j 32
    cd ..


## Setup Folder Structure (or reconfigure MS-EVAL/sv_util/settings.py)

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
    mkdir -p /MAdata/ena/simulated/UFRJ50816/Illumina-250
    mkdir -p /MAdata/ena/simulated/UFRJ50816/pacbio_CCS
    mkdir -p /MAdata/tmp


## Download Yeast Genomes & Build Indices

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


## Download Human Genome & Build Indices (will take a long time)

    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.27_GRCh38.p12/GCA_000001405.27_GRCh38.p12_genomic.fna.gz -O /MAdata/genome/human/GRCh38.p12/fasta/genome.fna.gz
    gunzip /MAdata/genome/human/GRCh38.p12/fasta/genome.fna.gz

    mkdir /MAdata/genome/human/GRCh38.p12/minimap
    ~/minimap2/minimap2 -x map-pb -d /MAdata/genome/human/GRCh38.p12/minimap/genome.map-pb.mmi /MAdata/genome/human/GRCh38.p12/fasta/genome.fna
    mkdir /MAdata/genome/human/GRCh38.p12/ngmlr
    cp /MAdata/genome/human/GRCh38.p12/fasta/genome.fna /MAdata/genome/human/GRCh38.p12/ngmlr/genome.fna
    mkdir /MAdata/genome/human/GRCh38.p12/ma
    ~/buildMA/maCMD --Create_Index /MAdata/genome/human/GRCh38.p12/fasta/genome.fna,/MAdata/genome/human/GRCh38.p12/ma,genome

## adjust RAM max usage of KSW

    adjust the parameter ksw_file_system_min_gb_size in MSV-Eval/sv_util/settings.py to match your ram - 4GB
## Run experiments

    cd MSV-EVAL
    python3 ambiguities_of_atomic_sv/main.py
    python3 svs_hidden_to_aligners/main.py
    python3 yeast/main.py



