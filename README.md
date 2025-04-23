
Person who wrote the code: Aline Muyle

contact: aline.muyle@cnrs.fr

Study: The contribution of small RNAs to the evolution of separate sexes and sex chromosomes in the plant Silene latifolia

Study summary:
Dioecy is a breeding system with separate females and males, where sex can be determined by sex chromosomes (for instance XY males and XX females). Dioecy is generally associated with the evolution of sexual dimorphism. In this study, we tested whether small RNAs (sRNAs) play a role in sexual dimorphism and sex chromosome evolution in Silene latifolia. We generated new data for female and male sRNAs and DNA methylation. We showed that sRNAs are most of the time female-biased in S. latifolia, suggesting the presence of the Y chromosome globally impacts the epigenome by diluting sRNAs genome-wide in males. We found limited evidence for the regulation of sex-biased genes by sRNAs, but we nonetheless identified a male-biased transcription factor that may potentially be regulated by sex-biased RNA-directed DNA methylation. This transcription factor might contribute to male traits, through the regulation of key factors in sex-determination and phenotypic sexual dimorphism. Finally, we compared female and male sRNA mapping along the S. latifolia sex chromosomes. We found that X and Y genes are targeted by significantly more sRNAs in males compared to females and PAR genes. Our results suggest that Y genes silencing following Y degeneration leads to the formation of sRNAs that can interact with both X and Y genes in males due to X-Y sequence homology. Our work calls for future investigation of the impact of these sRNAs generated from the Y chromosome on X gene expression in males. 

Code version 1.0

Overview of files and their contents:

-----------------------
epigenetics_XY_pairs.R
-----------------------
This R code allows to plot the graphs and run the GLMMs related to X/Y gene pairs. The input files listed below are necessary to run the code.

-------------
XY_pairs.txt
-------------
This is a table (tab delimited) listing X/Y gene pairs (one gene pair per line). 
The following columns can be found:
- geneX : the name of the X gene.
- chr : the chromosome of the gene.
- start : the start coordinate of the gene.
- end : the end coordinate of the gene.
- geneY : the name of the Y gene.
- strata : the strata of the gene pair.


-------------
genes.gff3
-------------
This is a simplified annotation file for the Silene latifolia genome v4.0 (Moraga et al. 2025), one line per gene.
The following columns can be found:
- chr : the chromosome of the gene.
- min : the start coordinate of the gene.
- max : the end coordinate of the gene.
- strand : the strand of the gene.
- GeneName : the name of the gene.


---------------------------
ptgs_sRNA_gene_mapping.tsv
---------------------------
This file lists sRNA mapping statistics for the PTGS pathway.
The following columns can be found:
- gene : the name of the gene.
- individual : the name of the individual.
- counts : the raw number of sRNA reads mapping.
- TPM : the normalized count of sRNA reads mapping in Transcripts per Million.
- tissue : the tissue sampled.
- sex : the sex of the individual.
- path : PTGS path (21 and 22nt sRNA).


---------------------------
rddm_sRNA_gene_mapping.tsv
---------------------------
Same as ptgs_sRNA_gene_mapping.tsv except for the path: rddm (24nt sRNA).


-------------------------------
y_v4_raw_gene_mapping_inds.tsv
-------------------------------
Raw numbers of sRNA reads mapping on Y genes (columns name similar to ptgs_sRNA_gene_mapping.tsv).


-------------------------------
y_v4_tpm_gene_mapping_inds.tsv
-------------------------------
Normalized counts of sRNA reads mapping (in TPM) on Y genes (columns name similar to ptgs_sRNA_gene_mapping.tsv).




