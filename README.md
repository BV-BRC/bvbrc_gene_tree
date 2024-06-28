# Gene Tree Service

## Overview

The BV-BRC Phylogenetic Tree Building Service enables construction of custom phylogenetic trees built from user-selected genomes, genes, or proteins. 
The standard interface makes it convenient to build a tree of a set of sequences (DNA or protein) that may come from a BV-BRC "feature group" (which can be built while exploring data on the website) or from user-uploaded fasta sequence files. 

The "FastTree" option computes large minimum evolution trees with profiles instead of a distance matrix. [1,2]. 
We also offer two maximum likelihood tree building algorithms: PhyML [3] and RaxML [4]. 
PhyML and RaxML infer a more evolutionarily accurate phylogenetic topology by applying a substitution model to the nucleotide sequences, but require the user to select an appropriate substitution model. 
This pipeline is best applied to datasets containing: 
1) fewer than 100 very long sequences, and
2) between 100 and 1,000 small or medium length sequences.

*Viral Genome Tree*
Another application of this pipeline is to small-to-medium non-segmented viral genomes.
This function is accessed through the "Viral Genome Tree" service interface.
Viral genomes of up to 250,000bp can be analyzed. The genomes are aligned by Mafft 

The service returns a Newick file which can be rendered in the interactive Archaeopteryx Tree Viewer in the BV-BRC or downloaded and viewed in other software.   



## About this module

This module is a component of the BV-BRC build system. It is designed to fit into the
`dev_container` infrastructure which manages development and production deployment of
the components of the BV-BRC. More documentation is available [here](https://github.com/BV-BRC/dev_container/tree/master/README.md).

This module provides the following application specfication(s):
* [GeneTree](app_specs/GeneTree.md)


## See also

* [Gene Tree Service Quick Reference](https://www.bv-brc.org/docs/quick_references/services/genetree.html)
* [Gene Tree Service](https://www.bv-brc.org/docs/https://bv-brc.org/app/GeneTree.html)
* [Gene Tree Tutorial](https://www.bv-brc.org/docs//tutorial/genetree/genetree.html)



## References

1.	Price MN, Dehal PS, Arkin AP. FastTree: computing large minimum evolution trees with profiles instead of a distance matrix. Mol Biol Evol. 2009 Jul;26(7):1641-50. doi: 10.1093/molbev/msp077. Epub 2009 Apr 17. PMID: 19377059; PMCID: PMC2693737. 
2.	Price MN, Dehal PS, Arkin AP. FastTree 2--approximately maximum-likelihood trees for large alignments. PLoS One. 2010 Mar 10;5(3):e9490. doi: 10.1371/journal.pone.0009490. PMID: 20224823; PMCID: PMC2835736.
3.	Guindon, S. and Gascuel, O., PhyML : "A simple, fast, and accurate algorithm to estimate large phylogenies by maximum likelihood." (2003) Syst Biol. 52: 696-704  
4.	Stamatakis, A. et al. (2005) Bioinformatics 21: 456-463
5.  Katoh K, Standley DM. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Mol Biol Evol. 2013 Apr;30(4):772-80. doi: 10.1093/molbev/mst010.
