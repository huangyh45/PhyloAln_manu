# PhyloAln_manu
## Data and scripts that reproduce the results and benchmark described in the manuscript "PhyloAln: a convenient reference-based tool to align sequences and high-throughput reads for phylogeny and evolution in the omic era"  
It should be noticed that the commands provided in each dataset folder are mainly used to provide the parameters for the analyses in the manuscript and may fail to directly rerun due to different file paths, installation and configuration. Detailed description of the procedures can be seen in the methods in the manuscript.  

### Directories
#### lifetree
Data and scripts used in the analyses of the simulated dataset across tree of life.  
#### contam
Data and scripts used in the analyses of the simulated dataset of contaminated fruit fly transcriptomes.  
#### ladybird
Data and scripts used in the analyses of the dataset of ladybird beetle transcriptomes.  
Raw data source: Li, et al. 2021. BMC Biology 19:7. http://doi.org/10.1186/s12915-020-00945-7  
#### plastome
Data and scripts used in the analyses of the dataset of pepper plastomes.  
Raw data source: Simmonds, et al. 2021. Molecular Phylogenetics and Evolution 163:107229. http://doi.org/10.1016/j.ympev.2021.107229  
#### UCE
Data and scripts used in the analyses of the dataset of turtle ultraconserved elements (UCEs).  
Raw data source: Crawford, et al. 2015. Molecular Phylogenetics and Evolution 83:250-257.  http://doi.org/10.1016/j.ympev.2014.10.021  

### Software requirements
Some tools are required in the commands:  
- PhyloAln (and its auxiliary scripts) v0.1.0 (https://github.com/huangyh45/PhyloAln/releases/tag/v0.1.0)
  
- gffread
- OrthoFinder
- MAFFT
- trimAl
- ART
- ReadSim
- Trinity
- CroCo
- BLAST
- EvidentialGene
- Read2Tree
- Orthograph
- IQ-TREE
- IDBA
