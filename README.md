# PhyloAln_manu
## Data and scripts that reproduce the results and benchmark described in the manuscript "PhyloAln: a convenient reference-based tool to align sequences and high-throughput reads for phylogeny and evolution in the omic era"  
It should be noticed that the commands provided in each dataset folder are mainly used to provide the parameters for the analyses in the manuscript and may fail to directly rerun due to different file paths, installation and configuration. Detailed description of the procedures can be seen in the methods in the manuscript.  

### File content
#### lifetree
Data and scripts used in the analyses of the simulated dataset across tree of life to conduct a comparative analysis of the alignment performance between PhyloAln and Read2Tree using different omic data, sequencing technologies, coverages, and species across tree of life.  
Files and directories:  
- CDS_trim: 46 single-copy gene alignments among 15 species using OrthoFinder + MAFFT + trimAl
- OGs_treelife_read2tree: Read2Tree result alignments
- phyloaln_nt_out: PhyloAln result alignments
- orthogroup.pl: (in ladybird folder) a custom script to run OrthoFinder
- orthomsa.pl: (in ladybird folder) a custom script to align the CDS sequences using OrthoFinder results
- readsim_multi.py: a custom script to run ReadSim in parallel
- run_commands.sh: the commands and parameters of the analyses in the dataset, can be opened in a text editor
#### contam
Data and scripts used in the analyses of the simulated dataset of contaminated fruit fly transcriptomes to evaluate capability of PhyloAln to eliminate foreign and cross cross-contamination and impact of contamination on Read2Tree.  
Files and directories:  
- CDS_trim: (in lifetree folder) 46 single-copy gene alignments among 15 species using OrthoFinder + MAFFT + trimAl
- OGs_fly_read2tree: Read2Tree result alignments
- phyloaln_b_add_nt_out: alignments of 'outgroup contamination' from the target species
- phyloaln_b_nt_out: PhyloAln result alignments without assembly step
- phyloaln_nt_out: PhyloAln result alignments with assembly step
- assemble.add_out.py: a modified assemble.py module to generate the alignments of 'outgroup contamination' from the target species
- assemble.debug: a modified assemble.py module to obtain the information of species source of the reads removed ot retained during decontamination
- contam.config: PhyloAln configure file
- phyloaln_b.aatree.rooted.tre: phylogeny of amino acid PhyloAln result alignments without assembly step
- phyloaln_b.speciestree.rooted.tre: phylogeny of CDS PhyloAln result alignments without assembly step
- phyloaln_b_add.aatree.rooted.tre: phylogeny of amino acid alignments of 'outgroup contamination' from the target species
- phyloaln_b_add.speciestree.rooted.tre: phylogeny of CDS alignments of 'outgroup contamination' from the target species
- run_commands.sh: the commands and parameters of the analyses in the dataset, can be opened in a text editor
- sample_reads_contam.py: a custom script to generate the simulated contamination dataset by randomly selected reads from the species sources
- sum_contam.py: a custom script to calculate the precision and recall of clean, foreign and cross contamination reads of each single-copy genes in PhyloAln alignments with assembly step
- sum_contam_b.py: a custom script to calculate the precision and recall of clean, foreign and cross contamination reads of each single-copy genes in PhyloAln alignments without assembly step
- sum_contam_outgroup.py: a custom script to calculate the average completeness and identity, and the precision and recall of clean, foreign and cross contamination reads of each single-copy genes in PhyloAln alignments with and without assembly step using different Drosophila species as the defined outgroups
#### ladybird
Data and scripts used in the analyses of the dataset of ladybird beetle transcriptomes to assess the performance of PhyloAln, Read2Tree and Orthograph on real transcriptome dataset.  
Raw data source: Li, et al. 2021. BMC Biology 19:7. http://doi.org/10.1186/s12915-020-00945-7  
Files and directories:  
- CDS_trim: 1290 single-copy gene alignments among 12 ladybird genomes
- orthograph_nt: Orthograph result alignments
- phyloaln_assembly_nt_out: PhyloAln result alignments using assemblies
- phyloaln_read_nt_out: PhyloAln result alignments using reads
- read2tree_result: Read2Tree result alignments
- single_copy_ref: 126 single-copy gene alignments among all 42 ladybird genomes and transcriptomes
- 00.config: PhyloAln configure file for 00 batch
- 01.config: PhyloAln configure file for 01 batch
- 02.config: PhyloAln configure file for 02 batch
- 03.config: PhyloAln configure file for 03 batch
- 04.config: PhyloAln configure file for 04 batch
- genome.list: the codes for 12 ladybird species with genomes
- orthograph.conf: Orthograph configure file, also required for rewriteconf.pl
- orthograph.og.tsv: OG member relationships for Orthograph
- orthograph.speciestree.rooted.tre: phylogeny of Orthograph result alignments
- orthogroup.pl: a custom script to run OrthoFinder
- orthomsa.pl: a custom script to align the CDS sequences using OrthoFinder results
- phyloaln_assembly.speciestree.rooted.tre: phylogeny of PhyloAln result alignments using the assemblies
- phyloaln_read.speciestree.rooted.tre: phylogeny of PhyloAln result alignments using the reads
- read2tree.speciestree.rooted.tre: phylogeny of Read2Tree result alignments
- rewriteconf.pl: a script to rewrite the Orthograph configure file and help loop to run Orthograph
- run_commands.sh: the commands and parameters of the analyses in the dataset, can be opened in a text editor
- seq.config: PhyloAln configure file for all the assemblies
- singlecopy.list: the list of 1290 single-copy genes among 12 ladybird genomes
- transassemble.pl: a custom script to assemble the transcriptomes  
##### Note: different batches represent the sequencing batches in which cross-contamination were detected and removed among the species by CroCo or PhyloAln
#### plastome
Data and scripts used in the analyses of the dataset of pepper plastomes to showcase the utility of PhyloAln in handling genes with non-standard genetic codes.  
Raw data source: Simmonds, et al. 2021. Molecular Phylogenetics and Evolution 163:107229. http://doi.org/10.1016/j.ympev.2021.107229  
Files and directories:  
- nt_out: PhyloAln result codon alignments
- test: de novo codon alignments using predicted genes
- plastid.config: PhyloAln configure file
- reftree.rooted.tre: phylogeny of de novo codon alignments using predicted genes
- run_commands.sh: the commands and parameters of the analyses in the dataset, can be opened in a text editor
- speciestree.rooted.tre: phylogeny of PhyloAln result alignments
#### UCE
Data and scripts used in the analyses of the dataset of turtle ultraconserved elements (UCEs) to evaluate the applicability of PhyloAln on nucleotide datasets such as UCEs.  
Raw data source: Crawford, et al. 2015. Molecular Phylogenetics and Evolution 83:250-257.  http://doi.org/10.1016/j.ympev.2014.10.021  
Files and directories:  
- all.fasta: the downloaded UCE matrix from Crawford, et al. (2015)
- aln.concatenated.fa: PhyloAln result alignment
- Graptemys_pseudogeographica.test.tsv: BLAST result of the assembly of Graptemys pseudogeographica reads to the sequences in the matrix
- reftree.rooted.tre: phylogeny of the downloaded UCE matrix
- run_commands.sh: the commands and parameters of the analyses in the dataset, can be opened in a text editor
- speciestree.rooted.tre: phylogeny of PhyloAln result alignment
- UCE.config: PhyloAln configure file
#### requirement.txt
The Conda configure file of major requirements.

### Requirements
#### Some tools are required in the commands:  
- gffread=0.12.1 (https://github.com/gpertea/gffread/releases/tag/v0.12.1)
- OrthoFinder=2.5.4 (https://github.com/davidemms/OrthoFinder/releases/tag/2.5.4)
- MAFFT=7.480 (https://mafft.cbrc.jp/alignment/software/source.html)
- trimAl=1.4.1 (https://github.com/inab/trimal/releases/tag/v1.4.1)
- ART=2016.06.05 (https://www.niehs.nih.gov/research/resources/software/biostatistics/art)
- Trinity=2.8.5 (https://github.com/trinityrnaseq/trinityrnaseq/releases/tag/Trinity-v2.8.5)
- BLAST=2.8.1 (https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.8.1/)
- Read2Tree=0.1.5 (https://github.com/DessimozLab/read2tree/releases/tag/v0.1.5)
- IQ-TREE=2.1.4_beta (http://www.iqtree.org/#download)
- IDBA=1.1.3 (https://github.com/loneknightpy/idba/releases/tag/1.1.3)  
The above tools can be installed through Conda and the Conda configure file provided in the repo, using the command:
```
conda install --file /your/PhyloAln_manu/path/requirement.txt
```
The rest tools can be manually installed:
- PhyloAln (and its auxiliary scripts) =0.1.0 (https://github.com/huangyh45/PhyloAln/releases/tag/v0.1.0)
- ReadSim=1.6 (https://sourceforge.net/projects/readsim/)
- CroCo=1.1 (https://gitlab.mbb.cnrs.fr/mbb/CroCo)
- EvidentialGene=2018.06.18 (http://arthropods.eugenes.org/EvidentialGene/other/evigene_old/evigene_older/)
- Orthograph=0.7.2 (https://github.com/mptrsen/Orthograph/releases/tag/0.7.2)
