# VisualizeRNAseq
R scripts used for visualization of RNA-seq data

These scripts were used for analysis of RNA-seq samples treated with a diverse collection of drugs.  To use these scripts, please carefully modify file paths to reflect the GSEA output of your local file structure, changing path names etc. 

## gatherGSEA
VisualizeRNAseq/gatherGSEA.R
The scripts assume a limited gene set across a large number of samples all in the same GSEA output folder.

### GSEA Bubble-Lattice Plot
The first visualization tool looks at Enrichment Scores for a set of gene sets vs. a set of samples:
![image](https://user-images.githubusercontent.com/11543307/41822725-56b817a6-77c2-11e8-906d-c3fea0ea44ac.png)




### GSEA Multi-sample Running Enrichment Plot
The second part of the script uses this same output but zooms in on a single gene set, looking at the shape of the running ES across multiple conditions:
![image](https://user-images.githubusercontent.com/11543307/41822705-1108769c-77c2-11e8-8238-a341f2dee81e.png)
