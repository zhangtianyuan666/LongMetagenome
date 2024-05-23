# Computational Tools and Resources for Metagenomics for Nanopore and PacBio sequencing

适用于长读长（PacBio和ONT）宏基因组软件和资源汇总  
更新时间(Update)：2024/5/23  
项目主页(Project homepage):https://github.com/zhangtianyuan666/LongMetagenome/  


## Figure 1
![image](https://github.com/zhangtianyuan666/LongMetagenome/assets/99855545/f10a21ec-11c0-4e00-8c06-3298e17aed41)  
Figure 1. The origin and development of the long-read metagenome study. Purple represents the origin stage of metagenomics, symbolizing the early beginnings and conception of the field (1991-2010). Green signifies the development of long-read metagenomics, indicating the period of advancement where longer DNA sequencing reads were introduced, enhancing the resolution and capabilities of metagenomic analysis (2011-2018).  Orange signifies the maturation and expansion phase of long-read metagenomics, highlighting a stage where this technology became more refined, widely adopted, and its applications broadened significantly (2019-2024+).   

## Figure 2 
![image](https://github.com/zhangtianyuan666/LongMetagenome/assets/99855545/2ecc11e2-bc75-4969-95ab-75c6abf080a0)  
Figure 2. Applying long-read metagenomics to analyze microbial community structure and functions.   

## Figure 3
![image](https://github.com/zhangtianyuan666/LongMetagenome/assets/99855545/a13f292c-83bb-4caa-afcc-4046fcf77352)  
Figure 3. Bioinformatics pipeline for analysis of long read metagenomic data.  

## Table 1 
<table>
  <tr>
    <th>Software</th>
    <th>Description</th>
    <th>Website</th>
  </tr>
  <tr>
    <td colspan="3">Data Quality Control, simulator, and remove host</td>
  </tr>
  <tr>
    <td>SMRTlink</td>
    <td>PacBio official workflows ranging from  base calling to sequence alignment et al</td>
    <td>https://www.pacb.com/support/software-downloads/</td>
  <tr>
    <td>bam2fastx</td>
    <td>Converts BAM formatted sequencing data into FASTQ format</td>
    <td>https://github.com/PacificBiosciences/bam2fastx</td>
  <tr>
    <td>Dorado</td>
    <td>A newer base-calling tool to instead Guppy</td>
    <td>https://github.com/nanoporetech/dorado</td>
  <tr>
    <td>PBSIM3</td>
    <td>A simulator for all types of PacBio and ONT long reads</td>
    <td>https://github.com/yukiteruono/pbsim3</td>
  <tr>
    <td>Porechop</td>
    <td>Adapter and chimera trimmer for Oxford Nanopore reads</td>
    <td>https://github.com/rrwick/Porechop</td>
  <tr>
    <td>NanoFilt</td>
    <td>Filtering and trimming of nonopore long-reads</td>
    <td>https://github.com/wdecoster/nanofilt</td>
  <tr>
    <td>LongQC</td>
    <td>Quality control of the PacBio and ONT long reads</td>
    <td>https://github.com/yfukasawa/LongQC</td>
   <tr>
    <td>Minimap2</td>
    <td>A versatile pairwise aligner for long-reads</td>
    <td>https://github.com/lh3/minimap2</td>
   <tr>
    <td>Winnowmap2</td>
    <td>Long-read or genome alignment software based on Minimap2</td>
    <td>https://gitlab.com/mcfrith/last</td>
   <tr>
    <td>LAST</td>
    <td>Pair-wise genome alignments</td>
    <td>https://gitlab.com/mcfrith/last</td>
   <tr>
    <td colspan="3">Taxonomy profile and binning</td>
  </tr>
    <td>Kraken2</td>
    <td>K-mer based taxonomic classifier</td>
    <td>https://ccb.jhu.edu/software/kraken2</td>
   <tr>
    <td>Bracken</td>
    <td>Bayesian estimation of abundance with Kraken</td>
    <td>https://ccb.jhu.edu/software/bracken/</td>
   <tr>
    <td>BugSeq</td>
    <td>Alignment and lower common ancestor (LCA) algorithm, highly accurate cloud platform for long-read metagenomic analyses</td>
    <td>https://bugseq.com/free</td>
   <tr>
    <td>Metamaps</td>
    <td>Mapping algorithm and expectation-maximization-based estimation for analysis of long-read metagenomic</td>
    <td>https://github.com/DiltheyLab/MetaMaps</td>
  <tr>
    <td>MEGAN-LR</td>
    <td>Alignment and LCA algorithm for taxonomic binning</td>
    <td>http://ab.inf.uni-tuebingen.de/software/downloads/megan-lr</td>
  <tr>
    <td>DESAMBA</td>
    <td>De Bruijn graph-based Sparse Approximate Match Block Analyzer (deSAMBA), a tailored long-read classification approach</td>
    <td>https://github.com/hitbc/deSAMBA</td>
  <tr>
    <td>Diamond</td>
    <td>Sequence aligner for protein and translated DNA searches, faster than BLAST</td>
    <td>https://github.com/bbuchfink/diamond</td>
   <tr>
    <td colspan="3">Metagenomic assemble, polish, and binning</td>
   <tr>
    <td>Hifiasm-meta</td>
    <td>Haplotype-resolved assembler for accurate HiFi reads</td>
    <td>https://github.com/lh3/hifiasm-meta</td>
   <tr>
    <td>metaFlye</td>
    <td>De novo assembler for long reads using repeat graphs</td>
    <td>https://github.com/fenderglass/Flye</td>
   <tr>
    <td>Lathe</td>
    <td>Generating bacterial genomes from metagenomes with Nanopore sequencing</td>
    <td>https://github.com/bhattlab/lathe</td>
   <tr>
    <td>metaMDBG</td>
    <td>Assembler for long and accurate metagenomics Reads (e.g. PacBio HiFi) based on the minimizer de-Brujin graph (MDBG)</td>
    <td>https://github.com/GaetanBenoitDev/metaMDBG</td>
 <tr>
    <td>STRONG</td>
    <td>Metagenomics strain resolution on assembly graphs</td>
    <td>https://github.com/chrisquince/STRONG</td>
 <tr>
    <td>Strainberry</td>
    <td>Automated strain separation of low-complexity metagenomes</td>
    <td>https://github.com/rvicedomini/strainberry</td>
 <tr>
    <td>OPERA-MS</td>
    <td>Hybrid metagenomic assembler which combines short and long-read</td>
    <td>https://github.com/CSB5/OPERA-MS</td>
 <tr>
    <td>Pilon</td>
    <td>Improves assemblies by correcting bases, fixing mis-assemblies and filling gaps via hierarchical polishing</td>
    <td>https://github.com/broadinstitute/pilon</td>
 <tr>
    <td>Racon</td>
    <td>Standalone consensus module to correct raw contigs via partial order alignment graph</td>
    <td>https://github.com/isovic/racon</td>
 <tr>
    <td>Medaka</td>
    <td>Correct draft sequences, create consensus sequences and variant calls from nanopore sequencing data via neural network model</td>
    <td>https://github.com/nanoporetech/medaka</td>
 <tr>
    <td>Ratatosk</td>
    <td>Hybrid error correction of long reads using colored de Bruijn graphs</td>
    <td>https://github.com/DecodeGenetics/Ratatosk</td>
 <tr>
    <td>MetaBAT2</td>
    <td>Similarity-based binner with label propagation algorithm</td>
    <td>https://bitbucket.org/berkeleylab/metabat</td>
 <tr>
    <td>metaWRAP</td>
    <td>Similarity-based binner with ensemble learning, integrated MetaBAT2, MaxBin2 and Concoct</td>
    <td>https://github.com/bxlab/metaWRAP</td>
 <tr>
    <td>metaBCC-LR</td>
    <td>Long reads binner with K-mer, composition, and density-based clustering</td>
    <td>https://github.com/anuradhawick/MetaBCC-LR</td>
 <tr>
    <td>LRBinner</td>
    <td>Long reads binner with K-mer and latent representation</td>
    <td>https://github.com/anuradhawick/LRBinner</td>
 <tr>
    <td>GraphMB</td>
    <td>Long reads binner with graph machine learning algorithms and the assembly graph generated</td>
    <td>https://github.com/MicrobialDarkMatter/GraphMB</td>
 <tr>
    <td>MetaCoAG</td>
    <td>Short and long reads binner via composition, coverage and assembly graphs</td>
    <td>https://github.com/metagentools/MetaCoAG</td>
 <tr>
    <td>MetaProb2</td>
    <td>Long reads binner with reads assembly, probabilistic k-mers, and graph modularity algorithm</td>
    <td>https://github.com/frankandreace/metaprob2</td>
 <tr>
    <td></td>
    <td>Row 3, Cell 2</td>
    <td>Row 3, Cell 3</td>
 <tr>
    <td>Row 3, Cell 1</td>
    <td>Row 3, Cell 2</td>
    <td>Row 3, Cell 3</td>
 <tr>
    <td>Row 3, Cell 1</td>
    <td>Row 3, Cell 2</td>
    <td>Row 3, Cell 3</td>
 <tr>
    <td>Row 3, Cell 1</td>
    <td>Row 3, Cell 2</td>
    <td>Row 3, Cell 3</td>
 <tr>
    <td>Row 3, Cell 1</td>
    <td>Row 3, Cell 2</td>
    <td>Row 3, Cell 3</td>
 <tr>
    <td>Row 3, Cell 1</td>
    <td>Row 3, Cell 2</td>
    <td>Row 3, Cell 3</td>
 <tr>
    <td>Row 3, Cell 1</td>
    <td>Row 3, Cell 2</td>
    <td>Row 3, Cell 3</td>
  </tr>
</table>
	


		
		

	Hybrid assembly and binning workflow for metagenomics, transcriptomics and pathway analysis	https://github.com/RVanDamme/MUFFIN

BASALT	Binning across a series of assembly’s toolkit for short and long read	https://github.com/EMBL-PKU/BASALT

HiCBin	Binning using Hi-C contact maps	https://github.com/dyxstat/HiCBin

MetaCC	Binning long-read and short-read metagenomic Hi-C data	https://github.com/dyxstat/MetaCC 

Nanodisco	Discovering multiple types of DNA methylation, and binning using nanopore sequencing	https://github.com/fanglab/nanodisco

dRep	Rapidly comparing large numbers of genomes and choosing the best representative genome	https://github.com/MrOlm/drep

GTDB-tk	Taxonomic classifications for bacterial and archaeal genomes	https://ecogenomics.github.io/GTDBTk/

Bugsplit	Highly accurate taxonomic binning of metagenomic assemblies	https://bugseq.com/academic

CheckM2	Predict the completeness and contamination of genomic bins using machine learning	https://github.com/chklovski/CheckM2

CoverM	Calculates coverage of genomes/MAGs	https://github.com/wwood/CoverM

metaQUAST	Evaluation of metagenome assemblies	http://bioinf.spbau.ru/metaquast

MetaCortex	Capturing variations in metagenomic assembly graphs	https://github.com/SR-Martin/metacortex

StrainPhlAn	Profiling microbes from known species with strain level resolution and providing comparative and phylogenetic	http://segatalab.cibio.unitn.it/tools/strainphlan/

MAGphase	Phasing for metagenomics using PacBio long reads	https://github.com/Magdoll/MagPhase

metaSVs	Combining long and short reads for analysis and visualization of structural variants in metagenomes	https://github.com/Wlab518/SV_procedure

Gene prediction and functional analysis	
Prokka	Rapid prokaryotic genome annotation	https://github.com/tseemann/prokka

HMMER	Searching sequence databases for sequence homologs by hidden Markov models	http://hmmer.org/

BLAST+	Basic Local Alignment Search Tool finds regions of similarity between biological sequences.	https://blast.ncbi.nlm.nih.gov/Blast.cgi

eggNOG-mapper	Functional annotation of novel sequences from the eggNOG database	http://eggnog-mapper.embl.de/

antiSMASH	Search a genome sequence for secondary metabolite biosynthetic gene clusters (BGCs)	https://antismash.secondarymetabolites.org/

BiG-SCAPE	Constructs sequence similarity networks of BGCs and groups them into cluster families	https://bigscape-corason.secondarymetabolites.org/

PlasFlow	Prediction of plasmid sequences in metagenomic contigs	https://github.com/smaegol/PlasFlow

PhiSpy	Finding prophages in bacterial genomes that combines similarity-and composition-based strategies	https://github.com/linsalrob/PhiSpy
Salmon	Highly-accurate, transcript-level quantification tools suitable metagenome	https://github.com/COMBINE-lab/salmon

Cd-hit	Clusters and compares protein or nucleotide sequences	https://github.com/weizhongli/cdhit

Table 1. The noteworthy software in long-read metagenomics studies.
