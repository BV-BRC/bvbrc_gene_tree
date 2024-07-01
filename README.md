# Gene Tree Service

## Overview

The BV-BRC Phylogenetic Tree Building Service enables construction of custom phylogenetic trees built from user-selected genomes, genes, or proteins. 

*Gene / Protein Tree*
The "Gene Tree" interface will build a tree of a set of gene sequences (DNA or protein) that may come from a BV-BRC "feature group" (which can be built while exploring data on the website) or from user-uploaded fasta sequence files (or a mixture of sources). 

*Viral Genome Tree*
The "Viral Genome Tree" interface will use the service to build a tree for entire viral genomes. 
It is limited to small-to-medium non-segmented viral genomes (up to 250,000 bp) of a single segment.

The sequences are aligned by MAFFT [1] unless a pre-aligned MSA file is submitted. 

The user has three choices for tree-building program: RAxML, PhyML, and FastTree.
RaxML [2] and PhyML [3] infer a rigorous maximum likelihood (ML) tree but take more time and require the user to select an appropriate substitution model. 
FastTree [4,5] computes minimum evolution trees with profiles in a way that approximates a ML tree but can handle much larger trees in shorter times.

The resulting tree is combined with user-specified metadata on each genome into a phyloxml file which can be viewed in the Archaeopteryx viewer [6] from the user's BV-BRC workspace.
Alternatively, the Newick version of the file can be downloaded and viewed in FigTree or other software.


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

1.  Katoh K, Standley DM. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Mol Biol Evol. 2013 Apr;30(4):772-80. doi: 10.1093/molbev/mst010.
2.	Stamatakis, A. et al. (2005) Bioinformatics 21: 456-463
3.	Guindon, S. and Gascuel, O., PhyML : "A simple, fast, and accurate algorithm to estimate large phylogenies by maximum likelihood." (2003) Syst Biol. 52: 696-704  
4.	Price MN, Dehal PS, Arkin AP. FastTree: computing large minimum evolution trees with profiles instead of a distance matrix. Mol Biol Evol. 2009 Jul;26(7):1641-50. doi: 10.1093/molbev/msp077. Epub 2009 Apr 17. PMID: 19377059; PMCID: PMC2693737. 
5.	Price MN, Dehal PS, Arkin AP. FastTree 2--approximately maximum-likelihood trees for large alignments. PLoS One. 2010 Mar 10;5(3):e9490. doi: 10.1371/journal.pone.0009490. PMID: 20224823; PMCID: PMC2835736.
6. Zmasek, Christian M.; Eddy, Sean R. (2001). "ATV: display and manipulation of annotated phylogenetic trees". Bioinformatics. 17 (4): 383–384. doi:10.1093/bioinformatics/17.4.383. PMID 11301314
