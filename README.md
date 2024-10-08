# Computational Tools and Resources for Long-read Metagenomic Sequencing using Nanopore and PacBio

适用于长读长（PacBio和ONT）宏基因组软件和资源汇总  
更新时间(Update)：2024/10/05
项目主页(Project homepage):https://github.com/zhangtianyuan666/LongMetagenome/  


## Figure 1
![Figure1_1](https://github.com/user-attachments/assets/6fa6641e-f25f-4cd6-8279-605f5546a907)
  Figure 1. The origin and development of the long-read metagenome studys. Purple represents the origin stage of metagenomics, symbolizing the early beginnings and conception of the field (1991-2010). Green signifies the development of long-read metagenomics, indicating the period of advancement where longer DNA sequencing reads were introduced, enhancing the resolution and capabilities of metagenomic analysis (2011-2018).  Orange signifies the maturation and expansion phase of long-read metagenomics, highlighting a stage where this technology became more refined, widely adopted, and its applications broadened significantly (2019-2024+).   

## Figure 2 
![Figure2_1](https://github.com/user-attachments/assets/bb704af8-e954-4cd5-a649-2bf890235b9c)
Figure 2. Applying long-read metagenomics to analyze microbial community structure and functions.   

## Figure 3
![Figure3_1](https://github.com/user-attachments/assets/13cd3b05-a691-498c-ba89-54077f2395b7)
Figure 3. Bioinformatics pipeline for analysis of long-read metagenomic data.  

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
    <td>A simulator for all types of PacBio and ONT long-reads</td>
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
    <td>Quality control of the PacBio and ONT long-reads</td>
    <td>https://github.com/yfukasawa/LongQC</td>
   <tr>
    <td>Minimap2</td>
    <td>A versatile pairwise aligner for long-reads</td>
    <td>https://github.com/lh3/minimap2</td>
   <tr>
    <td>Winnowmap2</td>
    <td>Long-read or genome alignment software based on Minimap2</td>
    <td>https://github.com/marbl/Winnowmap</td>
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
    <td>Alignment and lower common ancestor (LCA) algorithm, cloud platform for long-read metagenome</td>
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
    <td>De Bruijn graph-based Sparse Approximate Match Block Analyzer (deSAMBA), a tailored long-read classifier</td>
    <td>https://github.com/hitbc/deSAMBA</td>
  <tr>
    <td>Melon</td>
    <td>metagenomic long-read-based taxonomic identification and quantification using marker genes</td>
    <td>https://github.com/xinehc/melon</td>
  <tr>
    <td>Diamond</td>
    <td>Sequence aligner for protein and translated DNA searches, faster than BLAST</td>
    <td>https://github.com/bbuchfink/diamond</td>
 <tr>
    <td>metaBCC-LR</td>
    <td>Long-reads binner with K-mer, composition, and density-based clustering</td>
    <td>https://github.com/anuradhawick/MetaBCC-LR</td>
 <tr>
    <td>LRBinner</td>
    <td>Long-reads binner with K-mer and latent representation</td>
    <td>https://github.com/anuradhawick/LRBinner</td>
 <tr>
    <td>MetaProb2</td>
    <td>Long-reads binner with reads assembly, probabilistic k-mers, and graph modularity algorithm</td>
    <td>https://github.com/frankandreace/metaprob2</td>
   <tr>
    <td colspan="3">Metagenomic assemble, polish, and binning</td>
   <tr>
    <td>HiFiasm-meta</td>
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
    <td>Hybrid error correction of long-reads using colored de Bruijn graphs</td>
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
    <td>GraphMB</td>
    <td>Long-reads binner with graph machine learning algorithms and the assembly graph generated</td>
    <td>https://github.com/MicrobialDarkMatter/GraphMB</td>
 <tr>
    <td>MetaCoAG</td>
    <td>Short and long-reads binner via composition, coverage and assembly graphs</td>
    <td>https://github.com/metagentools/MetaCoAG</td>
<tr>
      <td>MUFFIN</td>
      <td>Hybrid assembly and binning workflow for metagenomics, transcriptomics and pathway analysis</td>
      <td>https://github.com/RVanDamme/MUFFIN</td>
    </tr>
    <tr>
      <td>BASALT</td>
      <td>Binning across a series of assembly’s toolkit for short and long-read</td>
      <td>https://github.com/EMBL-PKU/BASALT</td>
    </tr>
    <tr>
      <td>HiCBin</td>
      <td>Binning using Hi-C contact maps</td>
      <td>https://github.com/dyxstat/HiCBin</td>
    </tr>
    <tr>
      <td>MetaCC</td>
      <td>Binning long-read and short-read metagenomic Hi-C data</td>
      <td>https://github.com/dyxstat/MetaCC</td>
    </tr>
    <tr>
      <td>Nanodisco</td>
      <td>Discovering multiple types of DNA methylation, and binning using nanopore sequencing</td>
      <td>https://github.com/fanglab/nanodisco</td>
    </tr>
    <tr>
      <td>dRep</td>
      <td>Rapidly comparing large numbers of genomes and choosing the best representative genome</td>
      <td>https://github.com/MrOlm/drep</td>
    </tr>
    <tr>
      <td>GTDB-tk</td>
      <td>Taxonomic classifications for bacterial and archaeal genomes</td>
      <td>https://ecogenomics.github.io/GTDBTk/</td>
    </tr>
    <tr>
      <td>Bugsplit</td>
      <td>Highly accurate taxonomic binning of metagenomic assemblies</td>
      <td>https://bugseq.com/academic</td>
    </tr>
<tr>
      <td>CheckM2</td>
      <td>Predict the completeness and contamination of genomic bins using machine learning</td>
      <td>https://github.com/chklovski/CheckM2</td>
    </tr>
    <tr>
      <td>CoverM</td>
      <td>Calculates coverage of genomes/MAGs</td>
      <td>https://github.com/wwood/CoverM</td>
    </tr>
    <tr>
      <td>metaQUAST</td>
      <td>Evaluation of metagenome assemblies</td>
      <td>http://bioinf.spbau.ru/metaquast</td>
    </tr>
    <tr>
      <td>MetaCortex</td>
      <td>Capturing variations in metagenomic assembly graphs</td>
      <td>https://github.com/SR-Martin/metacortex</td>
    </tr>
    <tr>
      <td>StrainPhlAn</td>
      <td>Profiling microbes from known species with strain level resolution and providing comparative and phylogenetic</td>
      <td>http://segatalab.cibio.unitn.it/tools/strainphlan/</td>
    </tr>
    <tr>
      <td>MAGphase</td>
      <td>Phasing for metagenomics using PacBio long-reads</td>
      <td>https://github.com/Magdoll/MagPhase</td>
    </tr>
	 <tr>
    <td>Strainy</td>
    <td>Phasing and assembly of strain haplotypes using long-read data</td>
    <td>https://github.com/katerinakazantseva/strainy</td>
	</tr>
    <tr>
      <td>metaSVs</td>
      <td>Combining long and short reads for analysis and visualization of structural variants in metagenomes</td>
      <td>https://github.com/Wlab518/SV_procedure</td>
    </tr>
   <tr>
    <td colspan="3">Gene prediction and functional analysis</td>
   <tr>
    <tr>
      <td>Prokka</td>
      <td>Rapid prokaryotic genome annotation</td>
      <td>https://github.com/tseemann/prokka</td>
    </tr>
    <tr>
      <td>HMMER</td>
      <td>Searching sequence databases for sequence homologs by hidden Markov models</td>
      <td>http://hmmer.org/</td>
    </tr>
	<td>BLAST+</td>
      <td>Basic Local Alignment Search Tool finds regions of similarity between biological sequences</td>
      <td>https://blast.ncbi.nlm.nih.gov/Blast.cgi</td>
    </tr>
    <tr>
      <td>eggNOG-mapper</td>
      <td>Functional annotation of novel sequences from the eggNOG database</td>
      <td>http://eggnog-mapper.embl.de/</td>
    </tr>
    <tr>
      <td>antiSMASH</td>
      <td>Search a genome sequence for secondary metabolite biosynthetic gene clusters (BGCs)</td>
      <td>https://antismash.secondarymetabolites.org/</td>
    </tr>
    <tr>
      <td>BiG-SCAPE</td>
      <td>Constructs sequence similarity networks of BGCs and groups them into cluster families</td>
      <td>https://bigscape-corason.secondarymetabolites.org/</td>
    </tr>
    <tr>
      <td>PlasFlow</td>
      <td>Prediction of plasmid sequences in metagenomic contigs</td>
      <td>https://github.com/smaegol/PlasFlow</td>
    </tr>
    <tr>
      <td>PhiSpy</td>
      <td>Finding prophages in bacterial genomes that combines similarity-and composition-based strategies</td>
      <td>https://github.com/linsalrob/PhiSpy</td>
    </tr>
    <tr>
      <td>Salmon</td>
      <td>Highly-accurate, transcript-level quantification tools suitable for metagenome</td>
      <td>https://github.com/COMBINE-lab/salmon</td>
    </tr>
    <tr>
      <td>Cd-hit</td>
      <td>Clusters and compares protein or nucleotide sequences</td>
      <td>https://github.com/weizhongli/cdhit</td>
  </tr>
</table>
	

## Table 2
Table 2. The application of databases in long-read metagenomic studies

<table>
  <tr>
    <th>Database</th>
    <th>Description</th>
    <th>Tools</th>
    <th>Website</th>
  </tr>
  <tr>
    <td colspan="4">Functional annotation / reference databases</td>
  </tr>
  <tr>
    <td>Nr</td>
    <td>NCBI non-redundant database</td>
    <td>BLAST+</td>
    <td>https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/</td>
  </tr>
  <tr>
    <td>Uniprot</td>
    <td>Database of protein sequence and functional information for all species</td>
    <td>BLAST+</td>
    <td>https://www.uniprot.org/</td>
  </tr>
  <tr>
    <td>GO</td>
    <td>The Gene Ontology focuses on the function of the genes and gene products</td>
    <td>BLAST+, BLAST2GO</td>
    <td>https://www.geneontology.org/</td>
  </tr>
  <tr>
    <td>KEGG</td>
    <td>Kyoto Encyclopedia of Genes and Genomes</td>
    <td>Kofamscan, BLAST+, KOBAS</td>
    <td>https://www.genome.jp/kegg/</td>
  </tr>
  <tr>
    <td>Nt</td>
    <td>NCBI nucleotide database</td>
    <td>BLAST+</td>
    <td>https://www.ncbi.nlm.nih.gov/nucleotide/</td>
  </tr>
<td>RefSeq</td>
    <td>NCBI reference sequence database</td>
    <td>BLAST+</td>
    <td>https://www.ncbi.nlm.nih.gov/refseq/</td>
  </tr>
  <tr>
    <td>EggNOG</td>
    <td>Ortholog linkages, functional annotations, and gene evolutionary</td>
    <td>EggNOG-mapper</td>
    <td>http://eggnog5.embl.de/</td>
  </tr>
  <tr>
    <td>Rfam</td>
    <td>RNA families database</td>
    <td>HMMER</td>
    <td>https://rfam.org/</td>
  </tr>
  <tr>
    <td>Tigrfam</td>
    <td>Inferring protein families and domains based on Hidden Markov Models (HMMs)</td>
    <td>HMMER</td>
    <td>https://www.tigr.org/TIGRFAMs</td>
  </tr>
  <tr>
    <td>MBGD</td>
    <td>Microbial genome database for comparative analysis</td>
    <td>BLAST+</td>
    <td>https://mbgd.nibb.ac.jp/</td>
  </tr>
  <tr>
    <td colspan="4">Resistance and mobile genetic elements database</td>
  </tr>
  <tr>
    <td>mobileOG-db</td>
    <td>Bacterial mobile genetic elements</td>
    <td>BLAST+</td>
    <td>https://github.com/clb21565/mobileOG-db</td>
  </tr>
  <tr>
    <td>SARG 2.0</td>
    <td>Antibiotic resistance genes database</td>
    <td>ARGpore2, BLAST+, LAST</td>
    <td>http://smile.hku.hk/SARGs</td>
  </tr>
  <tr>
    <td>CARD</td>
    <td>The comprehensive antibiotic resistance database</td>
    <td>RGI, BLAST+</td>
    <td>https://card.mcmaster.ca/</td>
  </tr>
  <tr>
    <td>PHI</td>
    <td>Pathogen-host interactions</td>
    <td>BLAST+</td>
    <td>http://www.phi-base.org/</td>
  </tr>
  <td>VFDB</td>
    <td>Virulence factor database</td>
    <td>BLAST+</td>
    <td>http://www.mgc.ac.cn/VFs/</td>
  </tr>
  <tr>
    <td>BacMet</td>
    <td>Antibacterial biocide & metal resistance genes</td>
    <td>BLAST+</td>
    <td>http://bacmet.biomedicine.gu.se/</td>
  </tr>
  <tr>
    <td>MEGARes 3.0</td>
    <td>Resistance gene sequences for antimicrobial drugs, biocides, and metals</td>
    <td>BLAST+</td>
    <td>https://www.meglab.org/</td>
  </tr>
  <tr>
    <td>ISFinder</td>
    <td>Insertion sequences (IS) isolated from bacteria and archaea</td>
    <td>BLAST+</td>
    <td>https://isfinder.biotoul.fr/</td>
  </tr>
  <tr>
    <td>SecReT6 v3</td>
    <td>Type VI secretion system (T6SS)</td>
    <td>BLAST+</td>
    <td>https://bioinfo-mml.sjtu.edu.cn/SecReT6/</td>
  </tr>
  <tr>
    <td colspan="4">Metabolism and elemental cycling database</td>
  </tr>
  <tr>
    <td>CAZY</td>
    <td>Carbohydrate-active enZYmes database</td>
    <td>BLAST+, HMMER, dbCAN3</td>
    <td>http://www.cazy.org/</td>
  </tr>
  <tr>
    <td>CYPED</td>
    <td>Cytochrome P450 engineering database</td>
    <td>BLAST+</td>
    <td>http://www.cyped.uni-stuttgart.de</td>
  </tr>
  <tr>
    <td>TCDB</td>
    <td>Transporter classification (TC) system database</td>
    <td>BLAST+</td>
    <td>https://www.tcdb.org/</td>
  </tr>
  <tr>
    <td>antiSMASH</td>
    <td>Secondary metabolite biosynthetic gene clusters (BGCs)</td>
    <td>antiSMASH</td>
    <td>https://antismash.secondarymetabolites.org/</td>
  </tr>
  <tr>
    <td>Bigspace</td>
    <td>Diversity of biosynthetic gene clusters</td>
    <td>Bigspace</td>
    <td>https://bigscape-corason.secondarymetabolites.org/</td>
  </tr>
  <tr>
    <td>NCycDB</td>
    <td>Nitrogen cycle gene (sub) families</td>
    <td>BLAST+, Diamond</td>
    <td>https://github.com/qichao1984/Ncyc</td>
  </tr>
  <tr>
    <td>SCycDB</td>
    <td>Sulfur cycling gene and Pathways</td>
    <td>Diamond</td>
    <td>https://github.com/qichao1984/SCycDB</td>
  </tr>
  <td>MCycDB</td>
    <td>Methane cycling genes</td>
    <td>Diamond</td>
    <td>https://github.com/qichao1984/MCycDB</td>
  </tr>
  <tr>
    <td>PCyCDB</td>
    <td>Phosphorus cycling genes</td>
    <td>Diamond</td>
    <td>https://github.com/ZengJiaxiong/Phosphorus-cycling-database</td>
  </tr>
  <tr>
    <td>VB12Path</td>
    <td>Cobalamin synthesis pathways</td>
    <td>Diamond</td>
    <td>https://github.com/qichao1984/VB12Path</td>
  </tr>
  <tr>
    <td colspan="4">Taxonomic databases</td>
  </tr>
  <tr>
    <td>IMG/VR v4</td>
    <td>Integrated microbial genome/virus system</td>
    <td>BLAST+</td>
    <td>https://img.jgi.doe.gov/vr</td>
  </tr>
  <tr>
    <td>GTDB</td>
    <td>Genome taxonomy database</td>
    <td>GTDB-tk</td>
    <td>https://gtdb.ecogenomic.org/</td>
  </tr>
  <tr>
    <td>VirSorter2-DB</td>
    <td>Diverse DNA and RNA virus genomes</td>
    <td>VirSorter2</td>
    <td>https://github.com/jiarong/VirSorter2</td>
  </tr>
  <tr>
    <td>CheckV-DB</td>
    <td>Complete viral genomes from metagenomes</td>
    <td>CheckV</td>
    <td>https://bitbucket.org/berkeleylab/CheckV</td>
  </tr>
  <tr>
    <td>Kraken2-DB</td>
    <td>Standard or custom RefSeq databases for taxonomic classification</td>
    <td>Kraken2, Krakentools</td>
    <td>https://benlangmead.github.io/aws-indexes/k2</td>
  </tr>
  <tr>
    <td>Kaiju-DB</td>
    <td>Taxonomic classification database includes nr, RefSeq, progenomes, plasmid, and rvdb</td>
    <td>Kaiju</td>
    <td>https://bioinformatics-centre.github.io/kaiju/</td>
  </tr>
</table>


## Table 3
Table 3. The application of visualize R package in long-reads metagenomics studies  
| Application                             | Packages                                                                                                                                                        |
|-----------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Metagenome and microbiome analysis and visualization | [MetagenomeSeq](https://github.com/HCBravoLab/metagenomeSeq.git), [EasyAmplicon](https://github.com/YongxinLiu/Easyamplicon), [EasyMetagenome](https://github.com/YongxinLiu/EasyMetagenome), [EasyMicrobiome](https://github.com/YongxinLiu/EasyMicrobiome.git), [MicrobiomeStat](https://github.com/cafferychen777/MicrobiomeStat/tree/main), [microbiome](https://microbiome.github.io/tutorials/), [EasyMicroPlot](https://github.com/xielab2017/EasyMicroPlot.git), [Metacoder](https://grunwaldlab.github.io/metacoder_documentation/), [Phyloseq](https://joey711.github.io/phyloseq/)                                      |
| Data visualization & plotting           | [compositions](https://rdrr.io/cran/compositions/man/compositions.html), [bpca](https://www.rdocumentation.org/packages/bpca/versions/1.3-6/topics/bpca), [igraph](https://igraph.org/r/), [rtk](https://cran.r-project.org/web/packages/rtk/), [gbtools](https://github.com/Ensembl/gbtools), [Corrplot](https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html), [clusterProfiler](https://github.com/YuLab-SMU/clusterProfiler), [ImageGP](https://github.com/Tong-Chen/ImageGP)                                                                                     |
| Multi-omics                             | [ivTerm](https://github.com/SJTU-CGM/ivTerm/), [mixOmics](http://mixomics.org/)                                                                                                                                                 |
| Data processing & statistical analysis  | [ggplot2](https://www.rdocumentation.org/packages/ggplot2/versions/3.5.0), [ggthemes](https://www.rdocumentation.org/packages/ggthemes/versions/5.1.0), [ggtree](https://github.com/YuLab-SMU/ggtree), [ggmsa](https://github.com/YuLab-SMU/ggmsa), [networkD3](https://christophergandrud.github.io/networkD3/), [circlize](https://jokergoo.github.io/circlize/), [ggvenn](https://github.com/yanlinlin82/ggvenn), [ggmap](https://github.com/fresques/ggmap/), [ggpubr](https://rpkgs.datanovia.com/ggpubr/), [clusterCrit](https://github.com/cran/clusterCrit), [treemap](https://www.rdocumentation.org/packages/treemap/versions/2.4-4), [UpSetR](https://www.rdocumentation.org/packages/UpSetR/versions/1.4.0), [Pheatmap](https://www.rdocumentation.org/packages/pheatmap/versions/1.0.12/topics/pheatmap), [cowplot](https://github.com/wilkelab/cowplot)     |


# Supplementary Information
Table S1. The Software in Long-Read Metagenomics Studies.  
Table S2. The Description of Visualize R package in Long-Reads Metagenomics Studies.  
File S1. Installation and usage of the noteworthy metagenomic analysis software.  
