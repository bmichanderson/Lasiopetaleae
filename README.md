# Lasiopetaleae
Analyses of target sequence capture data for the tribe Lasiopetaleae (Malvaceae)  

This repository contains the scripts and steps for assembling and analysing Angiosperms353 and OzBaits data from the paper Anderson *et al.* 2025 "Target sequence capture informs generic delimitation and hybridization in the tribe Lasiopetaleae (Malvaceae)"  

Scripts are located in the `scripts` folder, an R markdown file in the `rmd` folder, and Singularity recipes in `singularity`  
Some additional label files and target files are in the `files` folder  

Raw sequencing data are available at the European Nucleotide Archive under projects PRJEB49212 (GAP stage 1) and PRJEB81152 (GAP stage 2)  
Note that uploading data to ENA removes Illumina information from read headers, so the filtering step to remove optical duplicates needs to be turned off  

The full analysis steps and notes can be followed in the `Lasio_analyses.md` file  
Many steps are specific to the NCI Gadi supercomputer environment  
