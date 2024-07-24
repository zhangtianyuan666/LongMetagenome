# File S1 Installation and usage of the noteworthy metagenomic analysis software

Update: 20240724

# Data Quality Control, simulator, and remove host
## SMRTlink `#f03c15`
Recommend using the official source code installation, it is easy to install

**Install**
```
# source code install
wget -c  https://downloads.pacbcloud.com/public/software/installers/smrtlink-release_13.1.0.221970.zip
unzip  smrtlink-release_13.1.0.221970.zip
./smrtlink-release_13.1.0.221970_linux_x86-64_libc-2.17_anydistro.run --rootdir smrtlink --smrttools-only
```
**Usage**
```
./smrtlink/smrtcmds/bin/ccs  subreads.bam  ccs.bam
```


## bam2fastx
Recommend using the conda to install

**install**
```
# conda
conad install bam2fastx
```
**Usage**
```
bam2fastq -o out.fq.gz -c 9 input.ccs.bam
```



## Dorado
Recommend using the official source code installation, it is easy to install

**Install**
```
#source code install
wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.6.2-linux-x64.tar.gz;
```
**Usage**
```
$dorado-0.5.3-linux-x64/bin/dorado basecaller -r --batchsize 3072 -x "cuda:all"  ./model/dna_r10.4.1_e8.2_400bps_sup@v4.2.0  ./pod5 > basecalled.bam
```


## PBSIM3
Recommend using the official GitHub installation, it is easy to install


**Install**
```
# 01. clone 
git clone https://github.com/yukiteruono/pbsim3.git
# 02. install
./configure;
make;
```
**Usage**
```
$pbsim3/src/pbsim --strategy wgs  --method sample  --sample sample/sample.fastq       --depth 20 --genome sample/sample.fasta       --sample-profile-id pf1
```



## Porechop
**Install**
The following three installation methods are acceptable,it is recommended to install using Conda first.
```
# conda install
conda install porechop
# docker 
docker pull quay.io/biocontainers/porechop:0.2.4--py36h2ad2d48_4
# singularity
singularity run https://depot.galaxyproject.org/singularity/porechop:0.2.4--py39hc16433a_4
```

Usage
```
porechop -i test.fq.gz \  # Input fastq
-o test.Adapter.fastq.gz \  # output fastq
-t 16 \ # threads
```

## NanoFilt
Recommend using the conda to install

**Install**
```
# conda
conda install -c bioconda nanofilt 
```
**Usage**
```
gunzip -c reads.fastq.gz | NanoFilt -q 10 | gzip > highQuality-reads.fastq.gz
```


## LongQC
Recommend using the conda and docker to install

**install**
```
# conda 
conda install longqc==1.2.0c
# docker
docker pull cymbopogon/longqc:1.2.0
# singularity
https://depot.galaxyproject.org/singularity/longqc%3A1.2.0c--hdfd78af_0
```
**Usage**
```
python longQC.py sampleqc -x ont-ligation -o out_dir input_reads.fq --ncpu 30
```

## Minimap2
Recommend using the conda to install

**Install**
```
# conda 
conda install minimap2
# docker 
docker pull nanozoo/minimap2:2.26--d9ef6b6
# singularity https://depot.galaxyproject.org/singularity/minimap2%3A2.9--1
```
**Usage**
```
# For PB data
minimap2 -ax map-ont -t 30 ref.fa   ont_read.fq > map.sam
# For ONT data
minimap2 -ax map-pb -t 30 ref.fa   pb_read.fq > map.sam
```

## Winnowmap2
Recommend using the conda to install

**Install**
```
# conda
conda install winnowmap==2.03
# docker
docker pull ctomlins/winnowmap:latest
# singularity
https://depot.galaxyproject.org/singularity/winnowmap%3A2.03--h5b5514e_1
```
**Usage**
```
meryl count k=15 output merylDB ref.fa
meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt
winnowmap -W repetitive_k15.txt -ax map-ont ref.fa ont.fq.gz > output.sam
```

## LAST
Recommend using the conda to install

**Install**
```
# conda 
 conda install last==1542
```
**Usage**
```
lastal -s 2 -T 0 -Q 0 -a 1 -P $N_threads -f BlastTab lastindex $Input_fa > blastresult
```


# Taxonomy Profile and binning
## Kraken2
Recommend using the conda to install

**Install**
```
# conda 
conda install kraken2==2.13
# docker
docker pull staphb/kraken2:2.1.3
# singularity
https://depot.galaxyproject.org/singularity/kraken2%3A2.1.3--pl5321hdcf5f25_0
```
**Usage**
```
kraken2 -db  minikraken2_v1_8GB/ -t 30 --report kraken2.report $fq --output kraken2.output
```

## Bracken
Recommend using the conda to install

**Install**
```
# conda 
conda install bracken
# docker
docker pull staphb/bracken:2.9
# singularity
https://depot.galaxyproject.org/singularity/bracken%3A2.9--py39h1f90b4d_0
```
**Usage**
```
bracken -d ${KRAKEN_DB} -i ${SAMPLE}.kreport -o ${SAMPLE}.bracken -r ${READ_LEN} -l ${LEVEL} -t ${THRESHOLD}
```
## Bugsplit/Bugseq
Recommend using online analysis

**Online analysis**
```
https://bugseq.com
```
**Install for local analysis**
```
# 01. install the nextflow
conda install nextflow
# 02. clone the packages
git clone https://github.com/awilson0/bugseq-pipeline.git
```
**Usage**
```
nextflow main.nf --fastq in.fq --outdir output_dir
```

## Metamaps
Recommend using the conda to install

**Install**
```
# conda
conda install metamaps
# docker
docker pull nanozoo/metamaps:latest
# singularity
https://depot.galaxyproject.org/singularity/metamaps%3A0.1.98102e9--h7ff8a90_1
```
**Prepare database**
```
# download a database, e.g. miniSeq+H and Extract the downloaded database into the databases/ directory
https://www.dropbox.com/s/6p9o42yufx2vkae/miniSeq%2BH.tar.gz?dl=0
```
**Usage**
```
# step1  map
metamaps mapDirectly -t 30 --all -r  ~benagen/Biodatabase/databases/miniSeq+H/DB.fa   -q SRR7497167_1.fastq.gz -o classification_results_metamaps_miniseq.SRR7497167

#step2 classify 
metamaps classify --DB /home/benagen/Biodatabase/databases/miniSeq+H/ --mappings classification_results_metamaps_miniseq.SRR7497167  -t 25 --minreads 10
```


## MEGAN-LR
Recommend using the Windows version

**Install**
```
# download the windows version
https://software-ab.informatik.uni-tuebingen.de/download/megan6/welcome.html
```

**Usage**
```
Either pipe the output of LAST directly to MAF2DAA, or apply MAF2DAA to the MAF file generated by LAST, to obtain a much smaller output file in DAA format.
Meganize the DAA file either using the Meganizer command-line tool or interactively in MEGAN.
```



## DESAMBA
Recommend using the official GitHub installation, it is easy to install

**Install**
```
# source code install
git clone https://github.com/hitbc/deSAMBA.git --depth=1
cd ./deSAMBA
bash ./build
bash ./build-index ./demo/viral-gs.fa ./demo_index
```

**Usage**
```
#classify
./bin/deSAMBA classify -t 4 ./demo_index ./demo/ERR1050068.fastq -o ./ERR1050068.sam
#analysis
./bin/deSAMBA analysis ana_meta_base ./ERR1050068.sam ./demo/nodes.dmp
```

## Diamond
Recommend using the official GitHub and conda installation , it is easy to install

**Install**
```
# conda 
conda install diamond
# Source code install
 http://github.com/bbuchfink/diamond/releases/download/v0.9.31/diamond-linux64.tar.gz # Source code install
tar xzf diamond-linux64.tar.gz # uncompress
# docker 
docker pull biocontainers/diamond-aligner:v0.9.24dfsg-1-deb_cv1
# singularity
https://depot.galaxyproject.org/singularity/diamond%3A2.1.9--h43eeafb_0
```
**Usage**
```
# build the database
diamond makedb --in test.fasta -d db_name
# align
diamond blastx -d db_name -q query.fasta -o output.txt
```

## metaBCC-LR
It is so difficult to install due to environmental conflicts!! Gcc version 9.4.0 needs to be installed on Ubuntu 20.04 LTS

**Install**
```
apt install libpthread-stubs0-dev
conda create -n metabcc python=3.8.19
 conda activate metabcc
pip install numpy scipy kneed seaborn h5py tabulate umap-learn song-vis
```
**Usage**
```
python mbcclr --resume -r test_data/data/reads.fasta -g test_data/data/ids.txt -o test_output -e umap -c 25000 -bs 10 -bc 10 -k 4
```

## LRBinner
Recommend using Docker/Singularity to run, followed by Conda installation

**Install**
```
# docker file
docker pull mantonius/lrbinner:v2.1
# source conda 
conda create -n lrbinner -y python=3.10 numpy scipy seaborn h5py hdbscan gcc openmp tqdm biopython fraggenescan hmmer tabulate pytorch pytorch-cuda=11.7 -c pytorch -c nvidia -c bioconda
# git clone 
git clone https://github.com/anuradhawick/LRBinner.git
# Build the binaries
sh build.sh
```
**Usage**
```
python lrbinners.py reads -r reads.fasta -bc 10 -bs 32 -o lrb --resume --cuda -mbs 5000 --ae-dims 4 --ae-epochs 200 -bit 0 -t 32
```
## MetaProb2
Recommend using Conda installation; MetaProb2 needs minimap2, Miniasm and MetaProb;

**Install**
```
# 01 Singularity
https://depot.galaxyproject.org/singularity/metaprob%3A2--boost1.61_1
# 02-1. Install minimap2, Miniasm
conda install minimap2;
conda install miniasm;
# 02-2. Source code install metaprob2
git clone https://github.com/frankandreace/metaprob2.git  #Download this repo
# 03 docker
docker pull flowcraft/metaprob:2-1
```
**Usage**
```
METAPROB2.sh [-s NUM SPECIES] [-k KMER-SIZE] [-w WINDOW-SIZE] [-m MAX-CHAINED-UTG-LENGTH] [-o OPT-PARAMETER-MODULARITY] [-l SKIP-READS-LEFT-OUT] [-r MIN-LENGTH MAX-LENGTH] <input_file> <output_folder> <name>
```



# Metagenome assemble, polish, and Binning
## Hifiasm-meta
Recommend using the official GitHub installation , it is easy to install

**Install**
```
# Source code install
git clone https://github.com/xfengnefx/hifiasm-meta.git
cd hifiasm-meta
make;
# docker image
docker pull bioinforpi/hifiasm-meta:v0.0.2
# singularity
https://depot.galaxyproject.org/singularity/hifiasm_meta%3Ahamtv0.3.1--h5b5514e_1
```
**Usage**
```
# assembly
hifiasm_meta -t32 -o asm reads.fq.gz 2>asm.log
hifiasm_meta -t32 -S -o asm reads.fq.gz 2>asm.log 
# Convert gfa to fasta
awk '/^S/{print ">"$2;print $3}' prefix.p_ctg.gfa > hifiasm.fa
```

## metaFlye
Recommend using  conda installation
**Install**
```
# conda install
conda install flye
# docker
docker pull quay.io/biocontainers/flye:2.9.3--py38he0f268d_1
# singularity
singularity run https://depot.galaxyproject.org/singularity/flye:2.9.3--py39hd65a603_1
```
**Usage**
```
flye --nano-raw no_Metazoa.pass.fq.gz --o
ut-dir flye_fq_20200101 --threads 40 --iterations 4 --meta
```



## Lathe
Recommend using conda installation

**Install**
```
# 01. install the snamemake
conda install snakemake
snakemake --version #please ensure this is >=5.4.3
# 02. git
git clone https://github.com/elimoss/lathe.git
```
**Usage**
snakemake --use-singularity --singularity-args '--bind /labs/,/scg/,/home/ ' -s /path/to/lathe/Snakefile \
--configfile path/to/modified_config.yaml --restart-times 0 --keep-going --latency-wait 30


## metaMDBG
Recommend using conda installation

**Install**
```
# conda
conda install -c conda-forge -c bioconda metamdbg
# docker
docker pull mmastert4/metamdbg:7ce646a
# singularity
https://depot.galaxyproject.org/singularity/metamdbg%3A0.3--hdcf5f25_0
```
**Usage**
```
metaMDBG asm ./path/to/assemblyDir reads.fastq.gz -t 4
```

## STRONG
Recommend using  conda installation

**Install**
```
./install_STRONG.sh
```
**prepare the gtdbtk database**
```
# Automatic: 
# 1. Run the command "download-db.sh" to automatically download and extract to: 
/opt/conda/envs/STRONG/share/gtdbtk-2.4.0/db/                      
# Manual: 
# 1. Manually download the latest reference data:
wget https://data.gtdb.ecogenomic.org/releases/release220/220.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r220_data.tar.gz         
# 2. Extract the archive to a target directory:
tar -xvzf gtdbtk_r220_data.tar.gz -C "/path/to/target/db" --strip 1 > /dev/null      rm gtdbtk_r220_data.tar.gz 
# 3. Set the GTDBTK_DATA_PATH environment variable by running:  
conda env config vars set GTDBTK_DATA_PATH="/path/to/target/db"
```
**Usage**
```
./bin/STRONG outputdir --config config.yaml
```


## Strainberry
Recommend using conda or docker images  

**Install**
```
# conda
git clone https://github.com/rvicedomini/strainberry.git
cd strainberry
conda env create -n sberry --file environment.yml
# docker
docker pull nanozoo/strainberry:v1.1--0e270f3
```
**Usage**
```
strainberry -r assembly.fasta -b alignment.sorted.bam -o sberry_out -c 4
```


## OPERA-MS
Recommend using docker images 

**Install**
```
# docker image
docker pull etheleon/opera-ms:latest
# Docker file
https://github.com/CSB5/OPERA-MS/blob/OPERA-MS-0.9.0/Dockerfile
#It integrated the following software: MEGAHIT,SPAdes,samtools,bwa,blasr,minimap2 ,Racon,Mash,MUMmer, and Pilon.
```
**Usage**
```
perl ../OPERA-MS.pl \ 
    --contig-file contigs.fasta \ #assembly contig by SR reads
    --short-read1 R1.fastq.gz \ #SR reads 
    --short-read2 R2.fastq.gz \ #SR reads 
    --long-read long_read.fastq \ #LR reads reads(Cannot use compressed gz format) 
    --no-ref-clustering \  #cluster
    --out-dir RESULTS \  #outdir
    2> log.err #log file
```
## Pilon
Recommend using conda installation

**Install**
```
# conda 
conda install pilon==1.24
# source code
git clone  https://github.com/broadinstitute/pilon/releases/download/v1.24/pilon-1.24.jar
# singularity
https://depot.galaxyproject.org/singularity/pilon%3A1.24--hdfd78af_0
# docker 
docker pull staphb/pilon:1.24
```
**Usage**
```
java -Xmn100G -Xms100G -Xmx250G -jar pilon-1.24.jar --genome Test.contigs.fasta --frags alignments.bam --fix all --outdir ./
```

## Racon
Recommend using Docker/Singularity to run, followed by Conda installation

**Install**
```
# conda
conda install racon
# docker
docker pull nanozoo/racon:1.5.0--c32e531
# singularity 
library://bioinfo-n1h/assemblatron/racon-1-5-0:latest
```
**Usage**
```
racon reads.fastq reads.paf assembly.fasta > corrected_assembly.fasta
# the paf file was producted by minimap2(map-pb or map-ont) 
minimap2 -ax map-pb assembly.fasta reads.fastq > reads.paf
```


## Medaka
Recommend using Docker/Singularity to run, followed by Conda installation

**Install**
```
# conda
conda install medaka
# docker
docker pull ontresearch/medaka:v1.11.3
# singularity
library://bioinfo-n1h/assemblatron/medaka-1-6-0:latest
```
**Usage**
```
medaka_consensus  -i NB01.spoa.fq -d NB01_trim.fastq -o NB01_spoa_medaka -t 4 -m  r1041_e82_400bps_hac_v4.2.0
```


## Ratatosk
We recommend using conda and Docker/Singularity installation and operation

**Install**
```
# conda
conda install ratatosk==0.9.0
# docker 
docker pull nanozoo/ratatosk:0.1--73685af
# singularity
https://depot.galaxyproject.org/singularity/ratatosk%3A0.9.0--hdcf5f25_0
```

**Uasge**
```
Ratatosk correct -v -c 30 -s M17_fastp.1.fq M17_fastp.2.fq -l M17_ont.fq -o M17_out_long_reads
```


## MetaBAT2
Recommend using Docker/Singularity to run, followed by Conda installation

**Install**
```
# conda install
conda install -c bioconda metabat2
# docker 
docker pull nanozoo/metabat2
# singularity
https://depot.galaxyproject.org/singularity/metabat2%3A2.15--h4da6f23_2
```
**Usage**
```
metabat2 -m 1500 \  # minimum_contig_length
-t 16 \  # threads
-i final.contigs.fa \  # contig_file
-a final.depth.txt \  # contig_depth
-o all \  #output_file
-v  # verbose_mode
```

## metaWRAP
Recommend using Docker/Singularity to run, followed by Conda installation

**Install**
```
# conda install
conda install -y -c bioconda metawrap
# docker 
docker pull quay.io/biocontainers/metawrap:1.2--0
# singularity
https://depot.galaxyproject.org/singularity/metawrap%3A1.2--hdfd78af_2
```
**Usage**
```
1.Bin the co-assembly
metawrap binning -o INITIAL_BINNING \
-t 96 \
-a ASSEMBLY/final_assembly.fasta \
--metabat2 \
--maxbin2 \
--concoct CLEAN_READS/ERR*fastq

2.Consolidate bin sets
metawrap bin_refinement -o BIN_REFINEMENT \
-t 96 \
-A INITIAL_BINNING/metabat2_bins/ \
-B INITIAL_BINNING/maxbin2_bins/ \
-C INITIAL_BINNING/concoct_bins/ \
-c 50 \
-x 10

3.Visualize
metawrap blobology -a ASSEMBLY/final_assembly.fasta \
-t 96 \
-o BLOBOLOGY \
--bins BIN_REFINEMENT/metawrap_50_10_bins \
CLEAN_READS/ERR*fastq

4.Re-assemble
metawrap reassemble_bins -o BIN_REASSEMBLY \
-1 CLEAN_READS/ALL_READS_1.fastq \
-2 CLEAN_READS/ALL_READS_2.fastq \
-t 96 \
-m 800 \
-c 50 \
-x 10 \
-b BIN_REFINEMENT/metawrap_50_10_bins

5.Classify_bins
metawrap classify_bins -b BIN_REASSEMBLY/reassembled_bins \
-o BIN_CLASSIFICATION \
-t 48

6.Annotate bins
metaWRAP annotate_bins -o FUNCT_ANNOT \
-t 96 \
-b BIN_REASSEMBLY/reassembled_bins/
```




## GraphMB
Recommend using Docker/Singularity to run
**Install**
```
# docker images
docker pull andrelamurias/graphmb:0.2.4 # images
# dockerfile
https://github.com/MicrobialDarkMatter/GraphMB/blob/main/Dockerfile 
# singularity
https://depot.galaxyproject.org/singularity/graphmb%3A0.2.5--pyh7cba7a3_0
```
**Usage**
```
graphmb --assembly data/strong100/ --outdir results/strong100/ --assembly_name edges.fasta --graph_file assembly_graph.gfa --depth edges_depth.txt --markers marker_gene_stats.tsv
```


## MetaCoAG
Recommend using Docker/Singularity to run, followed by Conda installation

**Install**
```
# conda
conda create -n metacoag -c bioconda metacoag
# docker 
docker://quay.io/biocontainers/metacoag:1.1.4--pyh7cba7a3_0
# singularity
https://depot.galaxyproject.org/singularity/metacoag%3A1.2.0--pyhdfd78af_0
```
**Usage**
```
metacoag --assembler flye --graph /path/to/assembly_graph.gfa --contigs /path/to/assembly.fasta --paths /path/to/assembly_info.txt --abundance /path/to/abundance.tsv --output /path/to/output_folder  
# assembler: spades,megahit,flye
```



## MUFFIN
Recommend using Conda installation, This is a nextflow pipeline.

**Install**
```
# 01.install the nextflow(Version 20.07+) by conda
conda install nextflow
# 02. git clone the pipeline
git clone https://github.com/RVanDamme/MUFFIN.git
```
**Usage**
```
# with conda
nextflow run RVanDamme/MUFFIN --output results_dir  --cpus 8 --memory 32g -profile local,conda,test
# with docker
nextflow run RVanDamme/MUFFIN --output results_dir  -profile local,docker,test
# MUFFIN comand:
nextflow run main.nf --ont nanopore/ --illumina illumina/ --assembler metaflye -profile conda
```



## BASALT
Recommend using Docker/Singularity to run, followed by Conda installation(difficult to install)

**Install**
```
# singularity
https://share.weiyun.com/xKmoBmrF

# official install
# clone
git clone https://github.com/EMBL-PKU/BASALT.git
cd BASALT
conda env create -n BASALT --file basalt_env.yml

# conda install
conda env create -n BASALT --file basalt_env.yml

# change permission to BASALT script
chmod -R 777 <PATH_TO_CONDA>/envs/BASALT/bin/*

#  Download the trained models for neural networks 'BASALT.zip' 
wget https://figshare.com/ndownloader/files/41093033
mv 41093033 BASALT.zip
mv BASALT.zip ~/.cache
cd ~/.cache
unzip BASALT.zip
```

**Usage**
```
# example1
BASALT -a as1.fa -s S1_R1.fq,S1_R2.fq/S2_R1.fq,S2_R2.fq -t 32 -m 128

# example2
BASALT\
-a as1.fa,as2.fa,as3.fa\
-s srs1_r1.fq,srs1_r2.fq/srs2_r1.fq,srs2_r2.fq\
-l lr1.fq,lr2.fq -hf hifi1.fq\
-t 60 -m 250\
--autopara sensitive --refinepara quick --min-cpn 40 --max-ctn 15 -qc checkm2
```

## HiCBin
We recommend using conda to run HiCBin.

**Install**
```
# 01. clone the repository 
git clone https://github.com/dyxstat/HiCBin.git
# 02. create a HiCBin environment using conda
cd HiCBin
conda env create -f HiCBin_linux_env.yaml
# 03. Download the R package
conda activate HiCBin_env
conda install r-glmmtmb

```
**Usage**
```
hicbin.py pipeline -e HpaII BN2WC.contigs.fa BN2WC1.sort.bam tax.txt BN2WC1.coverage.txt outdir
```


## MetaCC
We recommend using conda to install MetaCC

***Install**
```
# 01. Clone the repository
git clone https://github.com/dyxstat/MetaCC.git
# 02. Enter the MetaCC folder
cd MetaCC
# 03. Add dependencies of external softwares
chmod +x Auxiliary/test_getmarker.pl Auxiliary/FragGeneScan/FragGeneScan Auxiliary/FragGeneScan/run_FragGeneScan.pl Auxiliary/hmmer-3.3.2/bin/hmmsearch
# 04. Construct the conda environment
conda env create -f MetaCC_env.yaml
```
**Usage**
```
# step1 
MetaCC.py norm -e HpaII BN2WC.contigs.fa BN2WC1.sort.bam cover_dir
# step2
MetaCC.py bin --cover BN2WC.contigs.fa cover_dir
```

## Nanodisco
Recommend using Docker/Singularity to run

**Installation**
```
# Download the image from cloud.sylabs.io
singularity pull --name nanodisco.sif library://fanglab/default/nanodisco 
# Create a container named nd_env
singularity build nd_env nanodisco.sif          
```
**Prepare the container for examples**
```
# Create a writable container (directory) named nd_example
singularity build --sandbox nd_example nanodisco.sif 
# Start an interactive shell to use nanodisco, type `exit` to leave
singularity run --no-home -w nd_example 
```
**Usage**
```
# Methylation typing and fine mapping
# 01.get_data_bacteria # Retrieve E. coli current differences and reference genome
nanodisco characterize -p 4 -b Ecoli -d dataset/EC_difference.RDS -o analysis/Ecoli_motifs -m GATC,CCWGG,GCACNNNNNNGTT,AACNNNNNNGTGC -t nn -r reference/Ecoli_K12_MG1655_ATCC47076.fasta

# 02.Methylation binning of metagenomic contigs
nanodisco profile -p 4 -r reference/metagenome.fasta -d dataset/metagenome_subset_difference.RDS -w dataset/metagenome_WGA.cov -n dataset/metagenome_NAT.cov -b MGM1_motif -o analysis/binning --motifs_file dataset/list_de_novo_discovered_motifs.txt
nanodisco binning -r reference/metagenome.fasta -s dataset/methylation_profile_MGM1_motif.RDS -b MGM1_motif -o analysis/binning
nanodisco plot_binning -r reference/metagenome.fasta -u analysis/binning/methylation_binning_MGM1_motif.RDS -b MGM1_motif -o analysis/binning -a reference/motif_binning_annotation.RDS --MGEs_file dataset/list_MGE_contigs.txt
# Limit:  final_model_Albacore_v2.3.4_nn.RDS 
Albacore is not popular.
```


## dRep
The following three installation methods are acceptable,it is recommended to install using Conda first.

**Install**
```
# conda install
conda install drep
# docker
docker pull quay.io/biocontainers/drep:3.5.0--pyhdfd78af_0
# singularity
singularity run https://depot.galaxyproject.org/singularity/drep:3.5.0--pyhdfd78af_0
```
**Usage**
```
dRep dereplicate outout_directory -g path/to/genomes/*.fasta
```


## GTDB-tk
The following three installation methods are acceptable,it is recommended to install using Conda first.

**Install**
```
# conda install
conda install gtdbtk
# docker
docker pull quay.io/biocontainers/gtdbtk:2.4.0--pyhdfd78af_0
# singularity
singularity run https://depot.galaxyproject.org/singularity/gtdbtk:2.3.2--pyhdfd78af_0
```
**Usage**
```
gtdbtk  classify_wf --genome_dir ./ --out_dir classify_wf --extension fa --prefix bin --cpu 40
```

## Bugsplit/Bugseq
Recommend using online analysis

**Online analysis**
```
https://bugseq.com
```
**Install for local analysis**
```
# 01. install the nextflow
conda install nextflow
# 02. clone the git
git clone https://github.com/awilson0/bugseq-pipeline.git
```
**Usage**
```
nextflow main.nf --fastq in.fq --outdir output_dir
```


## CheckM2
Recommend using Docker/Singularity to run, followed by Conda installation

**Install**
```
# conda install
conda create -n checkm2 -c bioconda -c conda-forge checkm2
# docker 
docker pull marcoteix/checkm2:1.0.0
# singularity
https://depot.galaxyproject.org/singularity/checkm2%3A1.0.1--pyh7cba7a3_0
```
**Prepare database**
```
checkm2 database --download
```
**Usage**
```
checkm2 predict --threads 30 --input <folder_with_bins> --output-directory <output_folder> -x fa --threads  30
```


## CoverM
Recommend using Conda installation

**Install**
```
# conda install
conda install coverm
# singularity
https://depot.galaxyproject.org/singularity/coverm%3A0.7.0--h07ea13f_0
# docker 
docker pull biojokke/coverm:0.6.1
```
**Usage**
```
# 1.genome
coverm genome --coupled read1.fastq.gz read2.fastq.gz \  # input_fastq
--genome-fasta-files genome1.fna genome2.fna \  # genome_fasta_files
-o output.tsv  # output_file
# 2.contigs
coverm contig --coupled read1.fastq.gz read2.fastq.gz \  # input_fastq
--reference assembly.fna  # fasta_file_of_contigs
```

## metaQUAST
Recommend using Conda installation

**Install**
```
# conda install
conda install quast
# docker
docker pull cami/metaquast:latest
# singularity
https://depot.galaxyproject.org/singularity/quast%3A5.2.0--py39pl5321h4e691d4_3

```
**Usage**
```
metaquast.py -o quast -t 40 assembly.fasta
```

## MetaCortex
Recommend using Conda installation

**Install**
```
# conda
conda install -c bioconda metacortex
# singularity
https://depot.galaxyproject.org/singularity/metacortex%3A0.5.1--hec16e2b_1
```
**Usage**
```
echo "file1.fastq" > allreads.txt
echo "file2.fastq" >> allreads.txt
echo "file3.fastq" >> allreads.txt
metacortex_31 -k 31 -n 23 -b 65 -i allreads.txt -t fastq -f contigs.fa -l log.txt
```

## StrainPhlAn
Recommend using docker run

**Install**
```
# docker 
docker pull biobakery/strainphlan
```
**Usage**
```
strainphlan -s consensus_markers/*.pkl \  # SAMPLES
-m db_markers/t__SGB1877.fna \  # CLADE_MARKERS
-r reference_genomes/G000273725.fna.bz2 \  # reference_genomes
-o output \  # output_directory
-n 8 \  # threads 
-c t__SGB1877 \  # clade
--mutation_rates  # mutation_rates_table
```



## MAGphase
Recommend using the official GitHub installation, it is easy to install

**Install**
```
# clone the git 
git clone https://github.com/Magdoll/MagPhase.git
cd MagPhase
conda env create -f MagPhase.conda_env.yml
source activate MagPhase.env
python setup.py build
python setup.py install
```
**Usage** 
```
mag_phaser.py -a all_contigs.fasta -b align.sorted.bam -g gene.bed --bhFDR 0.01 -o 1377.strain
```


## metaSVs
Recommend using docker run

**Install**
```
# docker
docker pull wanglab518/metasvs:latest
```
**Prepare the database**
```
# 1.Downloading the CheckM database:
mkdir CheckM_db && cd CheckM_db
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
tar xvzf checkm_data_2015_01_16.tar.gz
rm checkm_data_2015_01_16.tar.gz

# 2.Downloading the gtdbtk database:
wget https://data.gtdb.ecogenomic.org/releases/release202/202.0/auxillary_files/gtdbtk_r202_data.tar.gz
tar xvzf gtdbtk_r202_data.tar.gz
mv release202 gtdbtk_db
rm gtdbtk_r202_data.tar.gz

# 3.Downloading the kobas database:
mkdir kobas_db && cd kobas_db
#to download sqlite3.tar.gz and seq_pep.tar.gz
Link: ftp://ftp.cbi.pku.edu.cn/pub/KOBAS_3.0_DOWNLOAD/
tar xvzf sqlite3.tar.gz  && tar xvzf seq_pep.tar.gz
rm sqlite3.tar.gz seq_pep.tar.gz
```
**Usage**
```
python /opt/conda/SV_procedure/call_SVs_procedure.py config.ini
```



# Gene prediction and functional analysis
## Prokka
This software is very easy to install, and it is recommended to use Conda. Additionally, Docker/Singularity can also be used for installation

**Install**
```
# conda install
conda install Prokka
# docker 
docker pull quay.io/biocontainers/prokka:1.14.6--pl5321hdfd78af_5
# singularity
singularity run https://depot.galaxyproject.org/singularity/prokka:1.14.6--pl5321hdfd78af_5
```
**Usage**
```
prokka --force --outdir output_prokka --prefix sample --metagenome --gcode 11 --cpus 10 --evalue 1e-5 --locustag locustag --addgenes --rnammer --rfam --compliant --centre sample contigs.fa
```


## HMMER
it is recommended to use Conda

**Install**
```
# conda install
conda install hmmer
# Source code installation
wget ftp://selab.janelia.org/pub/software/hmmer3/3.0/hmmer-3.0.tar.gz;
tar zxf hmmer-3.0.tar.gz;
cd hmmer-3.0;
    ./configure;
    make;
    make check;
# singularity
https://depot.galaxyproject.org/singularity/hmmer%3A3.4--hdbdd923_1
# docker
docker pull biocontainers/hmmer:v3.2.1dfsg-1-deb_cv1
```
**Usage**
```
# build the Hmm model
hmmbuild globins4.hmm tutorial/globins4.sto
# Using hmm search to search for protein sequence databases
hmmsearch globins4.hmm uniprot_sprot.fasta > globins4.out
```


## BLAST+
it is recommended to use Conda


**Install **
```
# conda
conda install blast

# source code install
adopt https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ Website Download
tar zxvpf ncbi-blast-2.10.1+-x64-linux.tar.gz 
 export PATH=$PATH:$HOME/ncbi-blast-2.10.1+/bin:
```
**Usage**
```
# build the database
makeblastdb -in xx.fasta -parse_seqids -hash_index -dbtype prot
# blastn
blastn -query seq.fasta -out seq.blast -db dbname -outfmt 6 -evalue 1e-5 -num_descriptions 10 -num_threads 8
# blastx
blastx -query seq.fasta -out seq.blast -db dbname -outfmt 6 -evalue 1e-5 -num_descriptions 10 -num_threads 8
```

## EGGNOG-mapper
This software is very easy to install, and it is recommended to use Conda. Additionally, Docker/Singularity can also be used for installation

**Install**
```
# conda install
conda install eggnog-mapper
# docker
docker pull quay.io/biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_0
# singularity
singularity run https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.12--pyhdfd78af_0
```
**Usage**
```
emapper.py -i $prokka_out/assembly.faa \ # cds
--cpu 40 \ # threads 
--data_dir $Public/Database/eggnog/  \ # database
--output assembly  \ # output
```

## AntiSMASH
This software is very easy to install, and it is recommended to use Conda. Additionally, Docker/Singularity can also be used for installation

**Install**
```
# conda
conda install antismash
# docker
docker pull quay.io/biocontainers/antismash:7.1.0--pyhdfd78af_0
# singularity
library://tamu_anand/default/antismash:latest
```
**Usage**
```
antismash --clusterblast --subclusterblast --knownclusterblast --smcogs --inclusive --borderpredict --full-hmmer --asf --enable-BiosynML --taxon bacteria --input-type nucl --cpus 20  sample1.fna
```

## BiG-SCAPE
This software is very easy to install, and it is recommended to use Conda. Additionally, Docker/Singularity can also be used for installation

**Install**
```
# conda
conda install bigscape
# docker 
docker pull currocam/big-scape:v1.1.5
# singularity
https://depot.galaxyproject.org/singularity/bigscape%3A1.1.9--pyhdfd78af_0
```
**Usage**
```
run_bigscape Streptomyces_BGC-MiBIG Streptomyces_BIGSCAPE_MiBIG -c 30 
```


## PlasFlow
This software is very easy to install, and it is recommended to use Conda. Additionally, Docker/Singularity can also be used for installation

**Inastall**
```
# conda install
conda install plasflow -c smaegol
# docker 
docker pull nanozoo/plasflow
# docker
https://depot.galaxyproject.org/singularity/plasflow%3A1.1.0--py35_0
```
**Usage**
```
PlasFlow.py --input Citrobacter_freundii_strain_CAV1321_scaffolds.fasta \
--output test.plasflow_predictions.tsv \
--threshold 0.7
```

## Phispy
This software is very easy to install, and it is recommended to use Conda. Additionally, Docker/Singularity can also be used for installation

**Install**
```
# conda
conda create -n phispy  -c bioconda phispy
# docker
docker pull nanozoo/phispy:3.7.8--d78a45e
# singularity 
https://depot.galaxyproject.org/singularity/phispy%3A4.2.21--py39h1f90b4d_3
```
**Usage**
```
PhiSpy.py genbank_file -o output_directory
```


## Salmon
The following three installation methods are acceptable,it is recommended to install using Conda first.

**Install**
```
# conda install
conda install salmon
# docker 
docker pull combinelab/salmon
# singularity
library://yh549848/rnaseqde/salmon:1.5.0
```
**Usage**
```
# 1.Build an index
salmon index -p 40 \  # Transcript fasta file
-t athal.fa.gz \  # salmon index
-i athal_index  # threads

# 2.Quant
salmon quant -i athal_index \  # salmon index
-l A  \  # type
-1 ${samp}_1.fastq.gz   \  #  input_fasta
-2 ${samp}_2.fastq.gz  \  #  input_fasta
-p 8  \  # threads 
--validateMappings \  # Quasi-mapping mode only
-o ${samp}_quant \ # output_file
--meta
```

## cd-hit 
This software is very easy to install, and it is recommended to use Conda. Additionally, Docker/Singularity can also be used for installation

**Install**

```
# conda
conda install cd-hit
# docker
docker pull biocontainers/cd-hit:v4.6.8-2-deb_cv1
# singularity
https://depot.galaxyproject.org/singularity/cd-hit%3A4.8.1--hdbcaa40_2
```
**Usage**
```
cd-hit-est -i All_sample.ffn \
-o cluster \
-c 0.95 \
-G 0 \
-aS 0.9 \
-g 1 \
-d 0 \
-T 32 \ #threads
-M 200000 #
```
