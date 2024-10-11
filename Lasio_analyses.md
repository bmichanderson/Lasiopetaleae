GAP Phylogenomics: Lasiopetaleae (Malvaceae)  
============================================
B.M. Anderson  
2023  

These notes describe the assembly and analysis of phylogenomic data for samples from the tribe Lasiopetaleae (Malvaceae)    

Scripts used are included in the folder `scripts` and include some scripts executed on a supercomputer (NCI Gadi) that contain path information and environmental variables not applicable to most users  

In some cases for commands run on Gadi, the shorthand `/path/to/...` is used in these notes and should be replaced with relevant paths to scripts and data being used  
It is assumed scripts are in a folder in the home directory `~/scripts`, or in a directory on the supercomputer `/path/to/scripts`  


# Data
The first stage of the Genomics for Australian Plants (GAP) Australian Angiosperm Tree of Life (AAToL) included the capture and sequencing of 9 samples from the group (sample numbers starting with 7 or 8)  
Bait-capture libraries were prepared using the Angiosperms353 bait kit (NEBNext Ultra II DNA with MyBaits Angiosperms353 Capture)  
Sequencing was done on an Illumina NovaSeq 6000 in the 150 bp paired-end read format  

For six samples (79869, 79870, 79871, 80849, 81208, 81210), an initial sequencing run was also done in the 100 bp paired-end read format  

For one sample (79740), library prep differed (LGC sbeadex oKtopure DNA extraction followed by NEBNext Ultra II DNA prep kit with Angiosperms353 MyBaits; the library was custom built with additional barcoding) and sequencing was on an Illumina NextSeq 500 in the 150 bp paired-end read format, trimmed to 143 bp  

As part of the second stage of AAToL, an additional 144 samples were sequenced (sample numbers starting with 3)  
Bait-capture libraries were prepared using the Angiosperms353 bait kit (NEBNext Ultra II FS DNA Library Prep with MyBaits Angiosperms353 Capture)  
Sequencing was done on an Illumina NovaSeq 6000 in the 150 bp paired-end read format  

To explore another bait kit, the same 144 samples from the second stage were captured and sequenced a second time  
Bait-capture libraries were prepared using the OzBaits bait kit (NEBNext Ultra II FS DNA Library Prep with MyBaits OzBaits Capture)  
Sequencing was done on an Illumina NovaSeq 6000 in the 150 bp paired-end read format  

## QC
Running FastQC on the raw data indicated issues in the first 9 bp of reads (or 1 bp in the pre-trimmed sample 79740), likely due to non-random enzymatic fragmentation  
A QC script (`illumina_qc.sh`) was run to remove the non-random bases (thought to potentially improve assembly), as well as remove optical duplicates, trim adapter sequences, and correct sequencing errors, keeping reads at least 50 bp long, using tools from the BBMap package  

To automate and parallelise the process, create a `calls.txt` file for each dataset  
The raw reads were named starting with the sample ID followed by an underscore (used for creating folders and naming QC output)  
Put a `samples.txt` file with sample IDs in each QC folder  
From the `/path/to/qc` folder with reads stored in a `/path/to/raw` folder:
```s
script="/path/to/scripts/illumina_qc.sh"
rawpath="/path/to/raw"
for sample in $(cat samples.txt); do
echo -e "${script} ${rawpath}/${sample}_*R1.fastq.gz ${rawpath}/${sample}_*R2.fastq.gz" >> calls.txt
done
```

Launch the jobs using a script to parallelise tasks
```s
qsub -l ncpus=192,mem=768GB,walltime=04:00:00,storage=gdata/nm31,wd -v calls_file="calls.txt",cores_per="16",timeout="4000" /path/to/scripts/launch_parallel.sh
```

Repeat for the stage1 samples in a folder `/path/to/qc/stage1` and the OzBaits in a folder `/path/to/qc/ozbaits`  

Following QC, stage1 reads of different original lengths (150 and 100 bp) were concatenated into single files for the six samples with two sets  


# Assembly
Multiple approaches to assembling the data were tried, including the SECAPR and HybPiper pipelines; the output from both could be formatted to be used as input to HybPhaser for evaluation and eventual phasing of putative hybrid samples  

For the data based on the Angiosperms353 bait set (most of the data), the first steps of the SECAPR pipeline could be run, but extracting target contigs had to be adapted to handle multi-exonic targets (multiple contigs assembled per target) using a custom script  
For the data based on the OzBaits bait set, the targets were single exons, so a single contig per target (in the absence of paralogs) was expected; in the case of multiple hits, stitching was run with the custom script  

## Target files
Two target files were created, one for each bait set (A353 and OzBaits)  

### Angiosperms353
To construct an initial set of targets for the Angiosperms353 data, the "NewTargets" (https://github.com/chrisjackson-pellicle/NewTargets) approach was used with their `mega353.fasta` file and the script `filter_mega353.py`, with the config file set to add "Malvaceae" targets; the resulting output was named `targets_Lasio.fasta`  

Since an annotated genome exists from a close relative (*Theobroma cacao*, GenBank accession GCF_000208745.1), CDS were extracted in protein format and used to construct a database to blast the targets and identify the associated CDS in *Theobroma*  
```s
# rename targets to simpler format
sed 's/ .*$//g' targets_Lasio.fasta | sed 's/_.*-/-/g' | sed 's/-.*-/-/g' > targets_renamed.fasta
# name the genbank flat file "Theo.gbff" and extract CDS as proteins
python3 ~/scripts/genbank_parse.py -l CDS Theo.gbff > CDS_list.txt
python3 ~/scripts/genbank_parse.py -f CDS_list.txt -t prot Theo.gbff
sed -i 's/ .*$//g' Theobroma_prot_extract.fasta		# remove extraneous text
python3 ~/scripts/remove_fastas.py "_" Theobroma_prot_extract.fasta		# remove extra copies (variants)
mv mod_Theobroma_prot_extract.fasta Theo_prot.fasta && rm Theobroma_prot_extract.fasta
# build database and blast targets (taking top hit)
makeblastdb -dbtype prot -in Theo_prot.fasta -out dbprot
blastx -query targets_renamed.fasta -db dbprot -num_threads 8 -outfmt 6 -max_target_seqs 1 > blast_out.tab
# sort results based on number of the first column (locus) and keep unique combos of locus and CDS
sed -E 's/^.*-(....\t)/\1/g' blast_out.tab | sort -k1 -n | awk '!seen[$1]++ && !seen[$2]++' | cut -f 1,2 > Theo_target_loci.tab
# extract the loci
cut -f 2 Theo_target_loci.tab > temp
python3 ~/scripts/genbank_parse.py -f temp -t nucl Theo.gbff
sed -i 's/ .*$//g' Theobroma_nucl_extract.fasta		# remove extraneous text
python3 ~/scripts/remove_fastas.py "_" Theobroma_nucl_extract.fasta		# remove extra copies (variants)
mv mod_Theobroma_nucl_extract.fasta Theo_nucl_A353.fasta && rm Theobroma_nucl_extract.fasta temp
# rename to the desired format (">taxon-locus") and use the A353 designation
while  IFS=$'\t' read -r -a myArray
do
locus353="${myArray[0]}"
locusTheo="${myArray[1]}"
sed -i "s/>${locusTheo}/>Theobroma-${locus353}/" Theo_nucl_A353.fasta
done < Theo_target_loci.tab
```
The resulting file `Theo_nucl_A353.fasta` contains only *Theobroma cacao* CDS sequences for 352 targets (locus 7602 was not found)  
After checking the target file with HybPiper, one locus (loc18592319 = locus 5354 in Angiosperms353) was detected as having internal stop codons; it had been adjusted in the ref seq but not when the locus was extracted; this was manually corrected by adding two "N"s at position 85 to correct the reading frame  

### OzBaits
The probe design can be traced back to exons from *Arabidopsis thaliana* genes  
Of the 100 target exons, 31 are from loci in the Angiosperms353 target set, so there is some overlap  

To create an initial target file, the specific exons were extracted from reference sequences for *Arabidopsis* chromosomes, then the main reading frame was used so that all targets were in-frame and without internal stop codons  
The NCBI accessions:
```
Chr1	CP002684.1  
Chr2	CP002685.1  
Chr3	CP002686.1  
Chr4	CP002687.1  
Chr5	CP002688.1  
```

Using a text file with the first field being the locus tag, the fifth being the chromosome number, and the sixth and seventh being the start and end positions (1-based and end-inclusive, with start > end indicating the other strand), the script `fasta_extract.py` can extract the exons  
Note: name the downloaded fasta files `Chr1.fasta`, `Chr2.fasta`, etc.  
```s
for locus in $(cut -f 1 exon_coords.tab); do
chrom="$(grep $locus exon_coords.tab | cut -f 5)"
coord="$(grep $locus exon_coords.tab | cut -f 6,7 | sed 's/\t/\.\./')"
python3 ~/scripts/fasta_extract.py -c $coord "Chr${chrom}.fasta"
sed -i "s/>.*/>Athaliana-$locus/" extract.fasta
cat extract.fasta >> exons.fasta
rm extract.fasta
done
```

The main reading frames were extracted using Biopython  
```python
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
with open('exons.fasta', 'r') as infile, open('exons_rf.fasta', 'w') as outfile1, open('exons_prot.fasta', 'w') as outfile2:
	exons = SeqIO.parse(infile, 'fasta')
	for exon in exons:
		longest_protseq = 0
		for frame in range(3):
			length = 3 * ((len(exon) - frame) // 3)
			nucleotides = exon.seq[frame: frame + length]
			running_coord = 0
			for protseq in nucleotides.translate().split('*'):
				if len(protseq) > longest_protseq:
					this_prot = protseq
					this_frame = frame
					this_nuc = nucleotides[running_coord: running_coord + 3 * len(protseq)]
					longest_protseq = len(protseq)
				running_coord = running_coord + 3 * len(protseq)
		outseq1 = SeqRecord(this_nuc, id = exon.id, description = exon.description, name = exon.name)
		outseq2 = SeqRecord(this_prot, id = exon.id, description = exon.description, name = exon.name)
		SeqIO.write(outseq1, outfile1, 'fasta')
		SeqIO.write(outseq2, outfile2, 'fasta')
```
The resulting `exons_rf.fasta` contains exon reading frames and `exons_prot.fasta` the peptide sequences for 100 exons  

Because SECAPR works better with closely-related sequences, similar hits were extracted from the *Theobroma cacao* genome  
Download the genome (as above for the A353 target file creation) and convert to fasta (`Theo.fasta`)  
Make a blast database, search with the protein sequences  
Keep hits > 50% identity  
```s
makeblastdb -dbtype nucl -in Theo.fasta -out dbnucl
tblastn -query exons_prot.fasta -db dbnucl -num_threads 8 -outfmt "6 qacc qlen sacc slen pident length qstart qend sstart send" -max_target_seqs 5 | awk '$5>50' > blast_out.tab
```
This reports hits for 98 of the 100 loci (missing AT2G30100 and AT5G04520)  

Retrieve the hits and write them with Biopython, in some cases stitching hits together in order  
```python
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
contigs = []
genome = SeqIO.parse(open('Theo.fasta', 'r'), 'fasta')
for contig in genome:
	contigs.append(contig)

blast_results = []
loci = []
with open('blast_out.tab', 'r') as blastfile:
	for line in blastfile:
		pieces = line.strip().split()
		locus = pieces[0].split('-')[-1]
		if locus not in loci:
			loci.append(locus)
		contig = pieces[2]
		pident = pieces[4]
		length = pieces[5]
		qstart = pieces[6]
		sstart = pieces[8]
		send = pieces[9]
		blast_results.append([locus, contig, pident, length, qstart, sstart, send])

loci = sorted(loci)
with open('Theo_nucl_ozbaits.fasta', 'w') as outfile, open('Theo_prot.fasta', 'w') as outfile2:
	for locus in loci:
		sublist = [entry for entry in blast_results if entry[0] == str(locus)]
		contig = sublist[0][1]
		subsublist = [entry for entry in sublist if entry[1] == contig]
		# sort by qstart and grab coordinates
		subsublist.sort(key = lambda x: int(x[4]))
		coords = [(int(entry[5]), int(entry[6])) for entry in subsublist]
		# create a fasta from the contig and coords
		myseq = Seq('')
		outname = 'Theobroma-' + str(locus)
		for contig_fasta in contigs:
			if contig_fasta.name == contig:
				this_contig = contig_fasta
				break
		for coord in coords:
			if coord[0] > coord[1]: 	# reverse strand
				myseq = myseq + this_contig.seq[coord[1] - 1: coord[0]].reverse_complement()
			else:
				myseq = myseq + this_contig.seq[coord[0] - 1: coord[1]]
		if len(myseq) < 100:
			print('No sequences long enough to retain for locus' + locus)
		else:
			new_record = SeqRecord(myseq, name = outname, id = outname, description = outname)
			prot_seq = Seq(str(myseq.translate()))
			new_record2 = SeqRecord(prot_seq, name = outname, id = outname, description = outname)
			SeqIO.write(new_record, outfile, 'fasta')
			SeqIO.write(new_record2, outfile2, 'fasta')
```
Two loci had sequences that were too short (AT2G38270 and AT5G01160)  
This leaves 96 loci in the target file `Theo_nucl_ozbaits.fasta`  

In one of the loci (AT1G15390), the stitched hits created an internal stop codon, so the codon was manually removed  

## SECAPR-like
SECAPR v. 2.2.8 was built in a Singularity container to run on Gadi  

The SPAdes assembly step of the SECAPR pipeline requires reads to be in a specified folder structure per sample  
Use symbolic links to the QC reads to create the required structure  
The sample numbers should be put in a text file (`samples.txt`), one per line  
```s
mkdir reads
for samplenum in $(grep ^3 samples.txt); do
mkdir reads/$samplenum
# for ozbaits, replace /path/to/qc with /path/to/qc/ozbaits
ln -s /path/to/qc/$samplenum/*.fastq.gz reads/$samplenum/
done
# for A353 only, with the stage 1 samples
for samplenum in $(grep ^[7,8] samples.txt); do
mkdir reads/$samplenum
ln -s /path/to/qc/stage1/$samplenum/*.fastq.gz reads/$samplenum/
done
```

Now launch the assembly script on Gadi that runs the equivalent step in the SECAPR pipeline  
```s
qsub -v qc_dir="reads" /path/to/scripts/secapr_assemble.sh
```
Total contigs averaged 27.9k per sample (5.6–165k; std dev 20.7k) for the A353 dataset, and 14.5k per sample (1.1–49.6k; std dev 9.3k) for the OzBaits data  
Note: the output folder `contigs/stats` contains non-essential assembly information and may be removed to save space and inode usage  

The target extraction step of the SECAPR pipeline requires a target file to use for blast searching  
Use the results of the previous step (folder `contigs`) and the target file (`Theo_nucl_A353.fasta` or `Theo_nucl_ozbaits.fasta`) as input  
```s
qsub -v contigs_dir="contigs",targets_file="Theo_nucl_A353.fasta" /path/to/scripts/secapr_targets.sh
```

For the Angiosperms353 data, the targets will have multiple matching contigs (flagged as "paralogs" by SECAPR), so selected contigs were stitched together separated by 100 "N"s by essentially following the logic:  
	Order the hits along the target and iterate through them  
	If a hit does not overlap with another hit (by more than 20 bp), keep the contig (if not already flagged for removal)  
	Else if two hits overlap, keep the one with the higher product of percentage identity and hit length (flag the other for removal)  

This was implemented in the script `secapr_stitch_contigs.py`  
To run, first link the extracted target contigs blast hits in a new directory  
The input to the script is a list of samples and the assembled contigs, and optionally the original names of the loci (SECAPR creates its own naming convention, so it needs to be converted back if wanting to keep track and consistent with HybPiper)  
```s
mkdir stitch_targets && cd stitch_targets
ln -s ../target_contigs/{3,7,8}*/{3,7,8}*_select* .
python3 /path/to/scripts/secapr_stitch_contigs.py -s ../samples.txt -t ../target_contigs/reference_fasta_header_info.txt ../contigs/*.fa
```
The output is a file for each sample named "sample"_targetcons.fasta in the `stitch_targets` folder  
The blast hits links can be removed (`rm stitch_targets/{3,7,8}*.txt`)  
Target contigs averaged 341 per sample (320–348; std dev 4.2)  

For the OzBaits data, the targets should ideally be hit by a single best contig in the absence of paralogs  
Target contigs were recovered for most loci, on average 93 per sample (92–95; std dev 0.7)  
Some "paralogs" were detected (contigs hitting the same target) (4–65 per sample; average of 20)  
In case some of the hits could be stitched, the contig stitching was run as with the A353 data  

## HybPiper
HybPiper v. 2.1.2 was run on Gadi  

For the Angiosperms353 data, an additional option (`--run_intronerate`) was added to recover supercontigs as well as exons  
For OzBaits data, HybPiper was run slightly differently (no `--run_intronerate` option used), as there was no need to account for introns  
Ensure there is a `samples.txt` file in the working directory with one sample per line  
Link reads and run the script (no `intronerate` option for OzBaits) with the respective target files for A353 or OzBaits  
```s
mkdir reads
for samplenum in $(grep ^3 samples.txt); do
# replace /path/to/qc with /path/qc/ozbaits for the OzBaits data
ln -s /path/to/qc/$samplenum/*.fastq.gz reads/
done
# only for A353
for samplenum in $(grep ^[7,8] samples.txt); do
ln -s /path/to/qc/stage1/$samplenum/*.fastq.gz reads/
done
# launch the job
qsub -v target_file="Theo_nucl_A353.fasta",reads_dir="reads",samples_file="samples.txt",target_type="dna",intronerate="y" /path/to/scripts/hybpiper.sh
```

The above jobscript will call HybPiper with the following arguments for each sample (with no `--run_intronerate` for OzBaits):
```s
hybpiper assemble --targetfile_dna Theo_nucl_A353.fasta -r {sample}_R1.fastq.gz {sample}_R2.fastq.gz --prefix {sample} --cpu 8 --diamond --merged --run_intronerate
```

Afterwards, to save file space, combine the results from all samples (launch from the directory with the `results` folder; use `intronerate="n"` for OzBaits)  
```s
qsub -v intronerate="y",delete="y" /path/to/scripts/combine_results.sh
```

For A353, samples had on average 348 loci with sequences (338–351; std dev 2.6), with an average of 274 loci (133–331; std dev 35.3) at 50% target length  
On average, 24% of reads were mapped (9.3–32.8; std dev 5.9)  
There were on average 0.7 (std dev 1.2) and 2 (std dev 3.9) loci per sample flagged as paralogs based on length and depth, with an outlier T. stelligera with 11 and 47  

For OzBaits, samples had on average 30 loci with contigs (24–35), with an average of 27 loci (23–32) at 50% target length  
On average, 1% of reads were mapped (0.4–1.5) suggesting problems with mapping to the targets  
There were on average 1.1 and 1.8 loci per sample flagged as paralogs based on length and depth, with an outlier Thomasia_stelligera with 18 and 18  

Since the recovery was poor for OzBaits, a second run was attempted  
To try to recover all the targets, HybPiper was run with additional arguments to make read mapping less stringent  
Arguments added: `--diamond_sensitivity sensitive` (> 40% identity; default: not set = > 60% identity) and `--evalue 0.001` (default: 0.0001)  
```s
qsub -v target_file="Theo_nucl_ozbaits.fasta",reads_dir="reads",samples_file="samples.txt",target_type="dna",add_args="--diamond_sensitivity sensitive --evalue 0.001" /path/to/scripts/hybpiper.sh
```

For the second OzBaits run, samples had on average 93 loci with sequences (92–94; std dev 0.6), with an average of 90 loci (87–92; std dev 0.7) at 50% target length  
On average, 18.5% of reads were mapped (8.1–27.4; std dev 4.3)  
There were on average 5 (std dev 8.9) and 8 (std dev 9.9) loci per sample flagged as paralogs based on length and depth, with an outlier T. stelligera with 74 and 77, and some other samples with high numbers (T. pauciflora with 53 and 55, L. ferraricollinum with 49 and 54)  


# Evaluation with HybPhaser
To assess the quality of the assemblies and the presence of hybrid samples and/or paralogous/problematic loci, HybPhaser v. 2.1 was run using the output supercontigs (which include "N"s between exons/contigs for proper read mapping), except for the OzBaits data assembled with HybPiper, which did not include the intronerate option  

HybPhaser was modified to address different file formatting and bug fixes in a fork of the GitHub repository (https://github.com/bmichanderson/HybPhaser/tree/develop_b); the modified version was built in a Singularity container to be run on Gadi  

While HybPiper output includes reads that were mapped for each locus, output from SECAPR target capture requires a mapping step to set up the files and folder structure that HybPhaser expects  

HybPhaser uses terminology as follows:  
SNPs = ambiguities in consensus sequences  
heterozygosity = percentage of loci with SNPs  
allelic divergence = proportion or percentage of positions in a locus that are SNPs, and analogously across all loci for a sample  

Recent hybrids from divergent taxa are likely to have high heterozygosity and high allelic divergence  
More ancient hybrid samples are unlikely to be detected with the evaluation step of HybPhaser, as they may not retain heterozygosity and divergent alleles  

To evaluate the data, flag loci recovered in < 75% of samples; for the remaining loci, flag outliers for excessive average allelic divergence (possible paralogs) with greater than 1.5 times the interquartile range above the third quartile (quartiles calculated based on only the loci recovered in >= 75% of samples)  

## SECAPR-like
The mapping job scripts can be launched with a parallel job script on Gadi and a `calls.txt` file that points to the inputs for each job  
The job script `hybphaser_readmap.sh` relies on a python script `split_gene_ref.py` to split combined output into individual loci files that HybPhaser expects; it runs mapping using tools from the BBMap package, bcftools and samtools  
```s
script="/path/to/scripts/hybphaser_readmap.sh"
reads_dir="/path/to/qc/"	# change this to /path/to/qc/ozbaits or /path/to/qc/stage1 as appropriate
for samplenum in $(cat samples.txt); do
echo -e "${script} ${reads_dir}${samplenum}/${samplenum}_R1.fastq.gz ${reads_dir}${samplenum}/${samplenum}_R2.fastq.gz stitch_targets/${samplenum}_targetcons.fasta" >> calls.txt
done
qsub -l ncpus=192,mem=384GB,jobfs=100GB,walltime=12:00:00,storage=gdata/nm31,wd -v calls_file="calls.txt",cores_per="8",timeout="7200" /path/to/scripts/launch_parallel.sh
```

The resulting output is a `HybPhaser` folder with the required input data (consensus sequences and mapped reads) to run evaluation and phasing steps  

The evaluation step of HybPhaser is run using a script that will automatically generate a `config.txt` file that can be modified for further analyses; as input, it requires the directory that was just produced  
The `intronerate` option needs to be selected given the way the files and folder structure have been named/generated  
No `seqlist` needs to be generated for the first evaluation step, as the data are not filtered yet  
(use `Theo_nucl_ozbaits.fasta` for the OzBaits dataset)  
```s
qsub -v out_dir="HybPhaser",target_file="Theo_nucl_A353.fasta",namelist="samples.txt",intronerate="y",seqlist="n" /path/to/scripts/hybphaser_evaluate.sh
```

For the A353 dataset:  
14 loci were recovered in < 75% of samples  
15 loci were average allelic divergence outliers, two of which were in the set recovered in < 75% of samples  
There was an evident break in the distribution of heterozygosity vs. allelic divergence; there were 23 samples with > 90% heterozygosity and > 0.6 allelic divergence (possible recent hybrids); one sample had 91.3% heterozygosity, but 0.46 allelic divergence and was on the non-hybrid side of the break (so was not included with the "hybrids")  

For the OzBaits dataset:
Three loci were recovered in < 75% of samples  
Four loci were average allelic divergence outliers  
There was an evident break in the distribution of heterozygosity vs. allelic divergence; there were 21 samples with > 85% heterozygosity and > 0.6 allelic divergence (possible hybrids)  

The two extra "hybrids" in the A353 dataset are stage 1 samples that weren't sequenced with the OzBaits baits  
Otherwise, the same 21 samples in both bait sets were detected as possible hybrids  

To collate the set of contig and consensus sequences for the desired loci (removing loci with low recovery and high allelic divergence), run HybPhaser and generate sequence lists after altering the `config.txt` file to point to a file of undesirable loci  
(change the parameter: `remove_loci_for_all_samples_with_more_than_this_mean_proportion_of_SNPs = "file"` and then set `file_with_putative_paralogs_to_remove_for_all_samples =` to path to `drop_loci.txt` which has one undesirable locus name per line)  
Change the name of the dataset in `config.txt` to `filtered`  
Keep the potential hybrids (filter them in later phylogenetic analysis)  
Set `remove_outlier_loci_for_each_sample= "yes"` in `config.txt` to also drop loci per sample when they are outliers for allelic divergence (potential gene duplications per sample rather than across samples)  
Do this for both data sets separately (with their respective target files, config files and folders)  
```s
qsub -v out_dir="HybPhaser",target_file="Theo_nucl_A353.fasta",namelist="samples.txt",intronerate="y",seqlist="y" /path/to/scripts/hybphaser_evaluate.sh
```

The `01_data` folder in the HybPhaser folder is large and contains many files, so after generating sequence lists, it can be tarred and compressed to save space and inode usage  

## HybPiper
The first step is to generate consensus sequences for evaluation of heterozygosity  
The HybPiper run produces a `results` folder, and the script running on Gadi compresses results for each sample, so a `decompress` option is needed if HybPiper was run with the script on Gadi  

When evaluating heterozygosity for bait sets with multiple exons per target, it is imperative that the HybPiper intronerated supercontigs are used (to avoid read mapping overlap across stitched exons); to make this more robust, the `hybphaser_generate.sh` script relies on the script `fasta_Ns_sub.py` that increases the number of "N"s used to connect contigs to 100 (HybPiper default is 10)  
The `bbmap` option ensures the consensus generation uses higher read mapping stringency than default (and in line with the SECAPR-like mapping step)  
(Use `intronerate="n"` for the OzBaits data)  
```s
qsub -v samples="samples.txt",hp2_dir="results",decompress="y",intronerate="y",bbmap="y" /path/to/scripts/hybphaser_generate.sh
```

As with the SECAPR-like approach, the next HybPhaser step uses the results of the first step for evaluation  
(use `Theo_nucl_ozbaits.fasta` and `intronerate="n"` for the OzBaits dataset)  
```s
qsub -v out_dir="HybPhaser",target_file="Theo_nucl_A353.fasta",namelist="samples.txt",intronerate="y",seqlist="n" /path/to/scripts/hybphaser_evaluate.sh
```

For the A353 dataset:  
Five loci were recovered in < 75% of samples  
19 loci were average allelic divergence outliers  
Note: three loci that were outliers in the SECAPR-like approach were not in these 19 (5791, 6198, 6128), and seven loci that were outliers here were not in the SECAPR-like approach (5944, 6977, 6785, 5802, 5816, 6407, 5404)  
There was an evident break in the distribution of heterozygosity vs. allelic divergence; there were 23 samples with > 95% heterozygosity (one sample with 93% but high allelic divergence, so a "hybrid") and > 0.9 allelic divergence (possible hybrids); one sample had 95.7% heterozygosity, but 0.64 allelic divergence and was on the non-hybrid side of the break (so was not included in the 23)  

These were the same 23 samples detected as "hybrids" in the SECAPR-like approach  
Overall, the allelic divergence was higher in the HybPiper assembly than in the SECAPR-like assembly  

For the OzBaits dataset:
Four loci were recovered in < 75% of samples  
Three loci were average allelic divergence outliers, all shared with the SECAPR-like approach  
Note: one locus that was an outlier in the SECAPR-like approach was not flagged as an outlier with HybPiper (AT1G75350)  
There was an evident break in the distribution of heterozygosity vs. allelic divergence; there were 21 samples with > 80% heterozygosity and > 0.9 allelic divergence (possible hybrids)  

These were the same 21 samples detected as "hybrids" in the SECAPR-like approach  
Allelic divergence was slightly higher in the HybPiper assembly compared to the SECAPR-like assembly  

Again, collate the set of contig and consensus sequences for the desired loci (removing loci with low recovery and high allelic divergence)  
Run HybPhaser and generate sequence lists after altering the `config.txt` file to point to a file of undesirable loci (`drop_loci.txt`) and changing the name of the dataset to `filtered`  
Set `remove_outlier_loci_for_each_sample= "yes"` in `config.txt` to also drop loci per sample when they are outliers for allelic divergence  
(use `Theo_nucl_ozbaits.fasta` and `intronerate="n"` for the OzBaits dataset)  
```s
qsub -v out_dir="HybPhaser",target_file="Theo_nucl_A353.fasta",namelist="samples.txt",intronerate="y",seqlist="y" /path/to/scripts/hybphaser_evaluate.sh
```

To determine which loci per sample were retained (for extraction of only those exons from HybPiper results rather than the intronerated contigs HybPhaser uses), tally from `HybPhaser/03_sequence_lists_filtered/samples_contigs`:  
```s
for file in *.fasta; do
grep ">" $file | sed 's/>//' | tr '-' '\t' >> sample_loci.tab
done
```
This produces a (long) file with two columns (sampleID and locus)  
Using this file and the results (`dna_seqs`) from HybPiper, files with exons from the filtered loci can be compiled in a folder `filtered_exons`  
```s
python3 /path/to/scripts/hybpiper_collect_loci.py -l ../sample_loci.tab /path/to/results/dna_seqs/*.fasta
```
This produced 328 files (exons only) for the A353 data and 89 files for the OzBaits data  


# Backbone trees
To infer backbone relationships and compare the two bait sets and assembly methods, phylogenetic analyses were run including all outgroups and excluding possible hybrid samples flagged during HybPhaser evaluation (23 samples for the A353 data; 21 of those also for the OzBaits data)  

For SECAPR-like data, input contigs for each locus were aligned with MAFFT v. 7.453 using the `--auto` mode, then cleaned with a script (`clean_alignment.py`) to remove positions with > 50% missing data and samples with > 75% missing data per locus  

For HybPiper data, input was restricted to only exons, so these were aligned using translated sequences before cleaning (as above) and phylogenetic analysis  

Tree inference was run in IQ-TREE v. 2.2.2, running a partitioned concatenation analysis of all loci (`-p {align_dir}`), and analyses of individual loci (`-S {align_dir}`), using a search for the optimal model (`-m MFP`), and 1000 Ultrafast Bootstrap replicates (`--ufboot 1000`), sampling by locus and site (`--sampling GENESITE`) in the concatenation analysis  
Branches in the locus trees with less than 50% bootstrap support were reduced to polytomies using Newick Utilities  
Gene and site concordance factors (using likelihood) were calculated using IQ-TREE and sampling 10,000 quartets for each internal branch for site concordance  
The loci trees were input to ASTRAL v. 5.7.1 to infer a tree under the multispecies coalescent, including a second run with a test for whether branches could be rejected as polytomies (`--branch-annotate 10`)  
The phylogenetic analyses were run on Gadi with the script `align_phylo.sh`  
The script relies on a `phylo.sif` Singularity container containing all the phylogenetic software  

In order to compare bait sets using identical sampling, stage 1 samples in the Angiosperms353 dataset were not included  

## SECAPR-like
For the A353 dataset, the input sequences were based on the filtered contigs from HybPhaser  
Samples had on average 315.4 (298–325) of 325 loci  
Dropping the 23 "hybrids" (includes two stage 1 samples) plus the remaining seven stage 1 samples (put all these in a `drop_samples.txt` file) leaves 123 samples  

For the OzBaits dataset, the input sequences were based on the filtered contigs from HybPhaser  
Samples had on average 87.7 (83–89) of 89 loci  
Dropping the 21 "hybrids" (put these in a `drop_samples.txt` file) leaves 123 samples  

For each dataset, phylogenetic analyses were run after excluding relevant samples  
From a folder called `phylo`:
```s
mkdir temp_input && cd temp_input
for locus_file in /path/to/HybPhaser/03_sequence_lists_filtered/loci_contigs/*.fasta; do
python3 /path/to/scripts/remove_fastas.py -f ../drop_samples.txt "$locus_file"
done
for file in mod*.fasta; do newname=$(echo $file | cut -f 2 -d "_"); mv $file ${newname}.fasta; done
cd ..
qsub -v align_dir="temp_input",realign="y",clean="y",analysis="a",poly="y" /path/to/scripts/align_phylo.sh
```

The filtered alignment for the A353 dataset consisted of 325 loci with a total length of 547,342 bp and 11.5% missing data  
The filtered alignment for the OzBaits dataset consisted of 89 loci with a total length of 70,706 bp and 9.4% missing data  

Make an `outgroup.txt` file with the sample numbers for all outgroup samples (*Seringia*, *Commersonia* and *Androcalva*)  
Create a `samples.tab` file with the sample IDs and sample names, tab separated, one per line  
Run scripts to convert the output trees to forms useful for plotting  
```s
python3 ~/scripts/concord_to_newick.py -t concord_gcf.cf.tree.nex -o concord_newick
python3 ~/scripts/concord_to_newick.py -t concord_scf.cf.tree.nex -o concord_newick
python3 ~/scripts/astral_parse.py -t astral.tre -f p -o astral
python3 ~/scripts/astral_parse.py -t astral.tre -f q -o astral
# If interested, can use Newick Utilities to collapse nodes with ASTRAL polytomy p-values > 0.05  
nw_ed astral_poly.tre "i & b > 0.05" o > astral_poly_collapsed.tre
```
Use the markdown file `plot_trees.rmd` interactively in R 4 for plotting  

## HybPiper
HybPiper output can be restricted to exons, potentially improving reliability of alignments when outgroups are included  
(This is not directly comparable to the SECAPR-like approach, where additional flanking regions may produce more usable characters and lead to greater tree resolution)  
Since the output is exons, these can be aligned using the translated sequences, then the nucleotide sequences can be aligned to correspond with the protein alignment  

For the A353 dataset, the input sequences were based on the filtered exons chosen after evaluation with HybPhaser  
Samples had on average 317.3 (302–328) of 328 loci  
Dropping the 23 "hybrids" (includes two stage 1 samples) plus the remaining seven stage 1 samples (put these in a `drop_samples.txt`) leaves 123 samples  

For the OzBaits dataset, the input sequences were based on the filtered exons chosen after evaluation with HybPhaser  
Samples had on average 88.3 (80–89) of 89 loci  
Dropping the 21 "hybrids" (put these in a `drop_samples.txt`) leaves 123 samples  

For a given set of exons, drop samples, translate to protein, substitute stop codons "*" with "X", align with MAFFT, then use `pal2nal.pl` (Suyama et al. 2006; http://www.bork.embl.de/pal2nal/distribution/pal2nal.v14.tar.gz) to convert the nucleotide files to align with the protein alignment  
(From a `phylo` folder where the `drop_samples.txt` file is located)  
```s
mkdir temp_input && cd temp_input
qsub -v dropfile="../drop_samples.txt",fasta_dir="/path/to/filtered_exons" /path/to/scripts/translate_align.sh
```

Remove the log files and any other unnecessary files before returning to the `phylo` folder  
Launch the phylogenetic analysis, cleaning alignments to remove positions with > 50% missing data and samples with > 75% missing data per locus (but no need for alignment this time)  
```s
qsub -v align_dir="temp_input",realign="n",clean="y",analysis="a",poly="y" /path/to/scripts/align_phylo.sh
```

The filtered alignment for the A353 dataset consisted of 328 loci with a total length of 255,159 bp and 7% missing data  
The filtered alignment for the OzBaits dataset consisted of 89 loci with a total length of 33,987 bp and 1.6% missing data  

As for the SECAPR-like approach, create relevant `outgroup.txt` and `samples.tab` files, modify the concordance and ASTRAL outputs for plotting, and use the markdown file `plot_trees.rmd` interactively in R 4 for plotting  

## Comparisons
As a concatenation and ASTRAL tree were generated for each combination of bait set (A353 and OzBaits) and assembly type (SECAPR-like and HybPiper), there were 2 * 2 * 2 = 8 tree comparisons  
Direct comparisons were made changing only one of either bait set or assembly type, fixing the other  

Assembler: SECAPR-like  
A353: 325 loci, 547 kbp  
OzBaits: 89 loci, 71 kbp  
A353 vs. OzBaits: congruent topologies for branching of major clades in ASTRAL trees, with one difference in concatenation trees (T. pygmaea sister to ingroup in A353 vs. main *Guichenotia* sister in OzBaits)  

Assembler: HybPiper (exons only)  
A353: 328 loci, 255 kbp  
OzBaits: 89 loci, 34 kbp  
A353 vs. OzBaits: largely congruent topologies for branching of major clades in ASTRAL trees, but more conflict in concatenation trees (T. stelligera on long branches and incongruent between bait sets; L. ferraricollinum on long branch in the OzBaits tree and placed sister to T. stelligera; Lys. abollatum separated from the other *Lysiosepalum* in the OzBaits tree)  

SECAPR-like vs. HybPiper:  
- for A353, the concatenation trees were largely congruent for the branching of major clades, while the ASTRAL trees differed in the placement of T. stelligera, with the SECAPR-like approach placing it as in the concatenation trees  
- for OzBaits, the ASTRAL trees were largely congruent except for the placement of T. stelligera, while the concatenation trees were less congruent within *Thomasia* and the placement of L. ferraricollinum and for the sister relationship of T. pygmaea or main *Guichenotia* to the rest of the ingroup (the HybPiper tree followed the A353 relationship)  

ASTRAL vs. concatenation: the ASTRAL trees tended to be more congruent for a given comparison  

Given the SECAPR-like approach produced more data compared to the way we extracted the HybPiper data, and given the problems associated with the T. stelligera sample in the HybPiper assemblies, we proceeded with the SECAPR-like approach  
Given the trees produced by the two bait sets are mostly congruent, we proceeded to combine the two bait sets (dropping loci in the OzBaits dataset that are shared with the A353 dataset)  
Further phylogenetic analyses need only focus on the broad ingroup (dropping *Commersonia*, *Androcalva* and *Seringia*), as it was recovered unambiguously as monophyletic in all the combinations of bait set and assembly/analysis approaches  


# HybPhaser for "hybrids"
To incorporate and phase the "hybrid" samples in further phylogenetic analyses of the ingroup, HybPhaser requires clade references to phase to  

Using the consistently recovered and well-supported clades in the backbone trees, select samples from them with low allelic divergence in the HybPhaser evaluation (least likely to be recent hybrids); references should also have good recovery, but most samples had good overall recovery of loci  
The choice of clade references is subjective and might require multiple iterations  
A challenging aspect of reference selection is ensuring divergent enough references to effectively phase reads unambiguously, while also not missing "hybrids" within more closely-related clades  
A second challenging aspect is that ancient hybrids or taxa that are on a shallow backbone (difficult to distinguish) may phase to multiple close clades because they share ancestral variation  

From the base of the ladderised tree, choose 14 references and set abbreviations:
- *Hannafordia*: 376791 H. bissillii subsp. latifolia (Hbis)
- *Thomasia pygmaea*: 376762 T. pygmaea (Tpyg)
- *Guichenotia* clade 1: 376672 G. tuberculata (Gtub)
- *Guichenotia* clade 2: 376668 G. micrantha (Gmic)

a divergent and odd *Guichenotia* clade:
- *Guichenotia* odd clade: 376658 G. anota (Gano)

the first Lasiopetalum clade (weakly supported):
- *Lasiopetalum* small clade: 376702 L. lineare (Llin)
- *Lasiopetalum* clade 1: 376799 L. micranthum (Lmic)
- *Lasiopetalum* clade 2: 376713 L. moullean (Lmou)

a second Lasiopetalum clade (strongly supported):
- *Lasiopetalum* clade 3: 376715 L. oldfieldii (Lold)  

a divergent *Thomasia* clade:
- *Thomasia julietiae* clade: 376754 T. julietiae (Tjul)

the main *Thomasia* clade (including *Lysiosepalum*):
- *Thomasia* small clade: 376779 T. sp. Hopetoun (TspH)
- *Lysiosepalum* clade: 376740 Lys. aromaticum (Laro)
- *Thomasia* remnant clade: 376749 T. foliosa (Tfol)
- *Thomasia* main clade: 376745 T. cognata (Tcog)

The sample numbers and abbreviations need to go into a csv file with the header `sample,abbreviation` called `clade_references.csv` in the `HybPhaser` folders  

## Clade association
The `config.txt` files need to be updated with a section specifying information for the clade association: a folder that will be created (`04_clade_association`), the csv file `clade_references.csv`, the reference consensus sequences folder, a folder to put the mapped reads and other settings depending on the computing resources  
Note that even though the original data is paired-end, the `single-end` setting should be used based on the way the reads have been collected (this is confusing but how the program works)  
Add these lines to the `config.txt` files, adjusting paths as appropriate:  
```s
# Clade association
path_to_clade_association_folder = "/path/to/HybPhaser/04_clade_association"
csv_file_with_clade_reference_names = "/path/to/HybPhaser/clade_references.csv"
path_to_reference_sequences = "/path/to/HybPhaser/03_sequence_lists_filtered/samples_consensus"
path_to_read_files_cladeassociation = "/path/to/HybPhaser/mapped_reads"
read_type_cladeassociation = "single-end"
ID_read_pair1 = ""
ID_read_pair2 = ""
file_with_samples_included = ""
path_to_bbmap = ""
no_of_threads_clade_association = "24"
run_clade_association_mapping_in_R = "no"
java_memory_usage_clade_association = "81G"
```
In addition, this step relies on the `01_data` folder (must be uncompressed/tarred)  

Launch the job from the folders where the `HybPhaser` folders are (with the updated `config.txt` files in them)  
```s
qsub -v out_dir="HybPhaser" /path/to/scripts/hybphaser_cladeassoc.sh
```

Afterwards, the folders `HybPhaser/mapped_reads` and `HybPhaser/ref` can be removed, and the `HybPhaser/01_data` folder is also no longer needed  

The output `04_clade_association/Table_clade_association.csv` gives the relative mappings to clade references for each sample  
While HybPhaser may identify more recent hybrids, some samples may represent ancient hybrids that possess a mixed genomic background with reduced heterozygosity but could still be impacting the resolution of the backbone of the tree  
To identify samples that are phasing to multiple major clades and may be impacting resolution of the backbone, open the clade associations in a spreadsheet program  
Combine sums of mapping for some close clade references that come from less resolved major clades (originally included to phase hybrids within groups):  
1) Thomasia group: TspH, Laro, Tfol, Tcog  
2) Lasiopetalum group: Llin, Lmic, Lmou  
3) Guichenotia group: Gtub, Gmic  

Then, for each sample, determine whether the top mapped reference (or reference group) comprises less than 60% of the total mapping, with at least 40% spread across other groups  
If so, that sample is potentially problematic for backbone support (possibly because of more ancient hybridisation), and might need to be excluded in a backbone tree  

Applying this threshold to the A353 results, nine ingroup samples were highlighted as potentially problematic, of which three are "hybrids"(*) between major groups:  
376674	L. angustifolium  
376673	L. adenotrichum  
376716	L. oppositifolium  
376688	L. ferraricollinum  
376681	L. compactum*  
376748	T. discolor*  
376736	L. sp. Desmond*  
80503	T. microphylla  
376772	T. stelligera  

Eight samples were highlighted for the OzBaits data, of which three are "hybrids"(*) between major groups:  
376673	L. adenotrichum  
376674	L. angustifolium  
376716	L. oppositifolium  
376681	L. compactum*  
376736	L. sp. Desmond*  
376688	L. ferraricollinum  
376748	T. discolor*  
376772	T. stelligera  

Some of these samples are on relatively isolated branches in the backbone trees with low support for placement, consistent with uncertainty caused by (ancient) admixture or lack of sampling of their relatives and diversity  
T. stelligera was already highlighted as potentially problematic during assembly with HybPiper, and it could represent contamination  
T. microphylla was only present in stage 1, so it had not yet been placed in a comparison tree  

To phase the samples that were detected as potentially recent hybrids by HybPhaser (see Evaluation above), the top clade references they associate with need to be identified for phasing  
Using the output file with the raw association information (`04_clade_association/Table_clade_association.csv`) along with the `clade_references.csv` file, run a simple script to choose the top references and prepare an output for the next phasing step in the proper format for HybPhaser  
In addition, select sample numbers excluding "hybrids" among the outgroup samples; this results in 16 of the 23 for the A353 data, or 15 of the 21 for the OzBaits data; put these in `hybrids.txt` files  
Choose to retain up to four references, as long as references had association values at least 33% of the top value  
From the `HybPhaser` folder:
```s
python3 /path/to/scripts/hybphaser_assess_references.py -c 04_clade_association/Table_clade_association.csv -r clade_references.csv -s hybrids.txt -p 0.33
```

If there are only two top references and both are from the same less-resolved clade, keep them (intra-clade "hybrid")  
If there are references from different major groups, keep them; if there are multiple references within less-resolved clades, only keep the top one when there are other referencs from different clades  
For example, if the top references for a sample are Llin and Lmou, keep them; if the top references are Llin, Lmou and Tjul, keep Llin and Tjul  

The association results are almost identical between the bait sets, but there are a couple of differences:  
376719	L. aff. parvuliflorum	A353 has Lmou second; OzBaits has Llin (Llin in A353 is almost the same, so pick Llin)  
376724	L. sp. RD 11032	A353 has Lmou second; OzBaits has Llin (Llin in A353 is almost the same, so pick Llin)  
376725	L. rosmarinifolium CW 1392	A353 has Lmou second; OzBaits has Llin (Llin in A353 is close, so pick Llin)  
376737	L. sp. Mt Ragged	A353 has Lmou second; OzBaits has Llin (Llin in A353 is almost the same, so pick Llin)  

376736	L. sp. Desmond	A353 has a third reference Tjul, while OzBaits only has the first two  
376748	T. discolor	OzBaits has a third clade reference Gano, while A353 only has two major clades; the top two references for both are Lmic and Lmou (intra-clade), but the order is different: choose Lmic for both  

In the two cases where one bait set has three top reference clades and the other has two, phase both bait sets to the shared two, but include a phasing to the ref in the other bait set even if it wasn't a top ref for that bait set  
Note that these two cases (with unique references between the two bait sets) may include non-biological processes such as contamination during library preparation and/or sequencing as an explanation for the difference  

Update the `phasing_prep.csv` file with manual adjustments as needed (e.g. dropping multiple references to the same clade, adding one where the other set has three)  

## Phasing
The phasing step requires the original read files, so create links in a `HybPhaser/read_links` folder  
```s
mkdir read_links
ln -s /path/to/qc/3*/*.gz read_links/		# OzBaits will have /path/to/qc/ozbaits/3*/*.gz  
ln -s /path/to/qc/stage1/{7,8}*/*.gz read_links/		# OzBaits will not have this
```

Update the `config.txt` file again to add a phasing section, specifying a folder to be created (`05_phasing`), the csv file just produced, the read links, the reference sequences to phase to, and where to put the phased reads (`phased_reads`)  
Note that now that the original read files are being used, the read type is `paired-end`  
Add these lines to `config.txt` (adjusting paths as necessary):  
```s
# Phasing
path_to_phasing_folder = "/path/to/HybPhaser/05_phasing"
csv_file_with_phasing_prep_info = "/path/to/HybPhaser/phasing_prep.csv"
path_to_read_files_phasing = "/path/to/HybPhaser/read_links"
read_type_4phasing = "paired-end"
ID_read_pair1 = "_R1.fastq.gz"
ID_read_pair2 = "_R2.fastq.gz"
reference_sequence_folder = "/path/to/HybPhaser/03_sequence_lists_filtered/samples_consensus"
folder_for_phased_reads = "/path/to/HybPhaser/phased_reads"
folder_for_phasing_stats = "/path/to/HybPhaser/05_phasing"
path_to_bbmap_executables = ""
no_of_threads_phasing = "24" 
java_memory_usage_phasing = "81G"
run_bash_script_in_R = "no"
```

Launch the jobs using the HybPhaser container from the directory where the `HybPhaser` folder is  
```s
qsub -v out_dir="HybPhaser" /path/to/scripts/hybphaser_phasing.sh
```

The folders `HybPhaser/ref` and `HybPhaser/read_links` can be removed afterward  
The `HybPhaser/05_phasing/Table_phasing_stats.csv` file provides a summary of the phasing to each reference  

## Assembly of phased reads
The phased reads now need to be assembled in a similar way to the original data (in this case, SECAPR-like)  
Naming convention needs to be adjusted and reads have to be formatted appropriately -- split into two -- for assembly; they also can be compressed to save space  

Run a simple job script on Gadi to format them appropriately (from the `HybPhaser` folder):
```s
qsub -v reads_dir="phased_reads" /path/to/scripts/hybphaser_phasedreadformat.sh
```

Set up the structure expected by SECAPR (folder for each sample) for the re-assembly  
```s
mkdir phased_links
for samplenum in $(ls phased_reads/ | cut -f 1 -d "_" | sort | uniq); do
mkdir phased_links/$samplenum
ln -s "$(pwd)"/phased_reads/${samplenum}_*.fastq.gz phased_links/$samplenum/
done
```

Launch the assembly script that runs the equivalent step in the SECAPR pipeline  
```s
qsub -v qc_dir="phased_links" /path/to/scripts/secapr_assemble.sh
```

Again, the resulting `contigs/stats` folder is unnecessary and can be removed to save space and inode usage  
Total contigs averaged about 1.2k per sample (835–2345) for the A353 dataset, and about 189 per sample (112–350) for the OzBaits data  

Collect the target contigs as done previously (use `../Theo_nucl_ozbaits.fasta` for OzBaits)  
```s
qsub -v contigs_dir="contigs",targets_file="../Theo_nucl_A353.fasta" /path/to/scripts/secapr_targets.sh
```

Average target recovery was ~310 targets per sample for the A353 data and ~87 targets per sample for the OzBaits data  

Stitch targets where appropriate  
```s
ls contigs/*.fa | sed 's/.*\///' | sed 's/\.fa//' > samples_phased.txt
mkdir stitch_targets && cd stitch_targets
ln -s ../target_contigs/3*/3*_select* .
ln -s ../target_contigs/7*/7*_select* .		# for A353 only
python3 /path/to/scripts/secapr_stitch_contigs.py -s ../samples_phased.txt -t ../target_contigs/reference_fasta_header_info.txt ../contigs/*.fa
```
The blast hits links can be removed (`rm {3,7}*.txt`) afterward

Note that moving forward, the assembled sequences (for phased "hybrids") need to be combined with the original data while dropping the corresponding original samples (unphased)  
To be consistent, loci that were previously dropped during HybPhaser filtering should be dropped again following phased assembly (if they were assembled)  

The `phased_links` and `phased_reads` folders are no longer needed  


# Combined phylogenetic analyses
The two bait sets can be combined to generate new trees, including a new "backbone" tree without questionable samples  
Now that the "hybrid" samples have been phased, they can also be included in a separate and more inclusive tree having all samples  
Ingroup samples (six) from stage 1 in the A353 dataset will not have the OzBaits loci (missing data)  

Use the SECAPR-like assembly approach for the combined analyses  
Make a new folder `combined_secapr` for collating the combined alignments and running analyses  

Overlapping loci in the OzBaits set already present in the A353 set (31) can be put in a `shared_loci.txt` file for dropping  

Of the 31 shared loci, 29 were recovered during OzBaits assembly, then three were dropped during HybPhaser evaluation of the OzBaits data, leaving 26 shared loci; if these are now excluded from the OzBaits data, that leaves 89 - 26 = 63 loci to add to the 325 from the A353 data, bringing the total to 388 loci  

To collate the sequences for the 388 loci, they need to be gathered from four sources:  
- A353 filtered loci  
- A353 phased loci ("hybrid" samples)  
- OzBaits filtered loci  
- OzBaits phased loci ("hybrid" samples)  

Create a `sequences` folder for collating  
Copy over the `samples.txt` file from the A353 dataset (all the samples)  
Create a `hybrids.txt` file with sample IDs of the 16 ingroup "hybrids", and an `outgroup.txt` file with sample IDs of 34 outgroup samples  

For the A353 filtered loci, link sequences keeping only the samples that are not "hybrids" or outgroups
```s
grep -vf <(cat hybrids.txt outgroup.txt) samples.txt > grab_samples.txt
for sample in $(cat grab_samples.txt); do
ln -s /path/to/A353_secapr/HybPhaser/03_sequence_lists_filtered/samples_contigs/"$sample"_* sequences/"$sample"_A353.fasta
done
```

For the A353 phased loci, the loci need to be filtered to keep only those in the corresponding original sample files (that had been filtered by HybPhaser)  
Make a file with names of all possible loci (352) using the original targets file: `grep ">" Theo_nucl_A353.fasta | cut -f 2 -d "-" > possible_loci.txt`  
```s
for sample in $(cat hybrids.txt); do
# link the hybrid phased fastas (multiple per hybrid)
for file in $(ls /path/to/A353_secapr/HybPhaser/stitch_targets/"$sample"*); do
filename="$(basename $file)"
ln -s "$file" sequences/"${filename/targetcons/A353_hyb}"
done
# determine what loci to drop, then drop them with a Python script
grep ">" /path/to/A353_secapr/HybPhaser/03_sequence_lists_filtered/samples_contigs/"$sample"_* | cut -f 2 -d "-" > keep_loci.txt
grep -v -f keep_loci.txt possible_loci.txt > drop_loci.txt
cd sequences
python3 /path/to/scripts/remove_fastas.py -f ../drop_loci.txt "$sample"*A353_hyb.fasta
cd ..
done
```

Rename the fastas without dropped loci to replace the links (from `sequences/`)  
```s
for file in mod_*; do
mv "$file" "${file/mod_/}"
done
```

Repeat this process for the OzBaits data **with adjusted paths**  
Make a file with names of all possible loci (96) using the original targets file: `grep ">" Theo_nucl_ozbaits.fasta | cut -f 2 -d "-" > possible_loci.txt`  
Changed commands from the above three blocks are (besides path updates):  
```s
...
grep -vf <(cat hybrids.txt outgroup.txt) <(grep ^3 samples.txt) > grab_samples.txt
...
for sample in $(grep ^3 hybrids.txt); do
...

```

As an additional step for the OzBaits data, remove shared loci (from `sequences/`)  
```s
python3 /path/to/scripts/remove_fastas.py -f ../shared_loci.txt *ozbaits*.fasta
for file in mod_*; do
mv "$file" "${file/mod_/}"
done
```

## New backbone
For the new backbone, drop the outgroups (*Seringia*, *Commersonia* and *Androcalva*), drop all "hybrid" samples, and drop the suspicious or conflicting samples that associated with multiple clade references (put IDs in `suspect_samples.txt`):  
376674	L. angustifolium  
376673	L. adenotrichum  
376716	L. oppositifolium  
376688	L. ferraricollinum  
80503	T. microphylla  
376772	T. stelligera 

Make a samples file for the backbone: `grep -vf <(cat hybrids.txt outgroup.txt suspect_samples.txt) samples.txt > backbone_samples.txt`  
This leaves 113 samples for the backbone tree  

Make a `phylo_backbone` folder, and a `temp_input` folder in that  
Run a script to collate the samples and loci, run from the `temp_input` folder 
```s
python3 /path/to/scripts/collate_loci_samples.py -s ../../backbone_samples.txt ../../sequences/*{ozbaits,A353}.fasta
```
This results in 388 locus files for alignment and analysis, each with up to 113 samples  

Run a first analysis from the `phylo_backbone` folder
```s 
qsub -v align_dir="temp_input",realign="y",clean="y",analysis="a",poly="y" /path/to/scripts/align_phylo.sh
```

The filtered alignment consisted of 388 loci with a total length of 604,731 bp and 11.3% missing data  

As for previous trees, create relevant `outgroup.txt` and `samples.tab` files, modify the concordance and ASTRAL outputs for plotting, and use the markdown file `plot_trees.rmd` interactively in R 4 for plotting  
In this case, the `outgroup.txt` should only have the three *Hannafordia* sample IDs  

It may be possible to improve the alignments and concordance factors by removing long branches (potential misassemblies) or by removing loci with little phylogenetic signal  

First, assess branch lengths for the loci using TreeShrink v. 1.3.9 (https://github.com/uym2/TreeShrink)  
Set the significance threshold to `0.10` rather than the default `0.05` to be slightly more sensitive at detecting outlier branches  
```s
singularity exec -H "$(pwd)" /path/to/phylo.sif run_treeshrink.py -t loci.treefile -m per-gene -q 0.10 -O output_ts -o treeshrink
```

Use the `temp_input` collated sequences (prior to alignment) but now drop TreeShrink outlier samples from the respective alignments  
```s
ls ../in_align/ | cut -f 1 -d "_" > loci.txt		# need to make sure the order is the same as in the first analysis
index=1
for locus in $(cat loci.txt); do
remove_line=$(sed -n "${index}p; $((index + 1))q" ../treeshrink/output_ts.txt | tr -s '\t' ',' | sed 's/,$//')
if [ ! -z "$remove_line" ]; then
python3 /path/to/scripts/remove_fastas.py "$remove_line" "${locus}.fasta"
else
echo "locus $locus does not have taxa to remove"
fi
((index+=1))
done
for file in mod_*; do
mv "$file" "${file/mod_/}"
done
```

Of the 113 samples, `cat ../treeshrink/output_ts.txt | tr -s '\t' '\n' | sort -n | uniq | wc -l` = 113 were dropped at least once  
Of the 388 loci, `grep -c "^$" ../treeshrink/output_ts.txt` = 110 had no drops  
Use `awk '{ print NF }' ../treeshrink/output_ts.txt` to tabulate drops per line; on average there were 2.9 drops per locus with a drop (std dev = 1.5; max = 8)  
All three of the outgroups were removed from six of the loci  

After removing samples from sequences, align and run phylogenetic analysis again, this time including likelihood mapping  
Make a folder for the run (`run2/`) and link the `temp_input` folder there  
```s
qsub -v align_dir="temp_input",realign="y",clean="y",analysis="a",poly="y",likemap="y" /path/to/scripts/align_phylo.sh
```

The filtered alignment consisted of 388 loci with a total length of 605,484 bp and 11.2% missing data  

Plot as before  

### Likelihood mapping
Assess the phylogenetic signal for loci using the results from likelihood mapping  

Likelihood mapping is a graphical tool to evaluate phylogenetic content of sequence alignments  
The method uses quartets of samples (three possible topologies) and assesses whether they are resolved or not by the data (i.e. is one topology favoured? Plot in one of the three triangle apices, areas 1,2, or 3)  
If all three topologies have roughly equal probability, then the quartet is not resolved and the resulting plotted point is in the centre region of the triangle (area 7)  
Area 7 represents "star"-like relationships vs. the more "treelike" regions (areas 1, 2 and 3)  
Likelihood mapping as implemented in IQ-TREE uses random choices (10,000 in our case) of four samples per alignment and assess how many of those quartets fall into the various areas of the triangle  

For each locus, a `*.lmap.quartetlh` file was produced along with log files reporting the results  
To count the numbers in area 7 (divide by 100 for percentage), use the `quartetlh` files  
Run from the `run2/` folder that has a `likemap` folder  
```s
ls in_align/ | cut -f 1 -d "_" > loci.txt
echo -e "Locus\tCount1\tCount2\tCount3\tCount4\tCount5\tCount6\tCount7" > likemap_counts.txt
for locus in $(cat loci.txt); do
counts=$(cut -f 8 likemap/lm_"$locus"*.lmap.quartetlh | sort -n | uniq -cd | tr -s ' ' | cut -f 2 -d ' ' | tr '\n' ' ')
paste <(echo $locus) <(echo $counts | tr ' ' '\t') >> likemap_counts.txt
done
```

Only 7 loci had > 40% of quartets in area 7, suggesting most loci were phylogenetically informative  
On average, loci had 17.4% of quartets in area 7 (std dev = 9.1; max = 59.0; min = 0.8)  
No loci had > 11% of quartets in the partly resolved areas (4, 5 and 6)  

Quartets in area 7 may indicate lack of resolving power for parts of the tree, but the locus may still be informative for other parts of the tree  
There doesn't seem to be a good reason to drop loci from the analysis under suspicion they are lowering concordance because of phylogenetic noise or lack of resolution  

## All samples
For the inclusive tree, drop the outgroups (*Seringia*, *Commersonia* and *Androcalva*), but include the phased copies of the "hybrid" samples as well as the suspect samples  

Make a samples file, adding the new phased sample IDs (34 of them; e.g. `79740Tcog` for the `79740` sample phased to reference abbreviated `Tcog`); put these IDs in a `phased_samples.txt` file  
```s
grep -vf <(cat hybrids.txt outgroup.txt) samples.txt > temp
cat temp phased_samples.txt > master_samples.txt && rm temp
```
This leaves 153 samples for analysis  

Make a `phylo_all` folder, and a `temp_input` folder in that  
Run a script to collate the samples and loci run from the `temp_input` folder  
```s
python3 /path/to/scripts/collate_loci_samples.py -s ../../master_samples.txt ../../sequences/*.fasta
```
This results in 388 loci files for alignment and analysis, each with up to 153 samples  

Run a first analysis from the `phylo_all` folder
```s 
qsub -v align_dir="temp_input",realign="y",clean="y",analysis="a",poly="y" /path/to/scripts/align_phylo.sh
```

The filtered alignment consisted of 388 loci with a total length of 581,675 bp and 12.5% missing data  

As for previous trees, create relevant `outgroup.txt` and `samples.tab` files, modify the concordance and ASTRAL outputs for plotting, and use the markdown file `plot_trees.rmd` interactively in R 4 for plotting  
Again, the `outgroup.txt` should only have the three *Hannafordia* sample IDs  

As with the backbone tree, investigate removing long branches (potential misassemblies) and loci with little phylogenetic signal  

Assess branch lengths for the loci using TreeShrink  
```s
singularity exec -H "$(pwd)" /path/to/phylo.sif run_treeshrink.py -t loci.treefile -m per-gene -q 0.10 -O output_ts -o treeshrink
```

Use the `temp_input` collated sequences (prior to alignment) and drop TreeShrink outlier samples from the respective alignments  
```s
ls ../in_align/ | cut -f 1 -d "_" > loci.txt		# need to make sure the order is the same as in the first analysis
index=1
for locus in $(cat loci.txt); do
remove_line=$(sed -n "${index}p; $((index + 1))q" ../treeshrink/output_ts.txt | tr -s '\t' ',' | sed 's/,$//')
if [ ! -z "$remove_line" ]; then
python3 /path/to/scripts/remove_fastas.py "$remove_line" "${locus}.fasta"
else
echo "locus $locus does not have taxa to remove"
fi
((index+=1))
done
for file in mod_*; do
mv "$file" "${file/mod_/}"
done
```

Of the 153 samples, `cat ../treeshrink/output_ts.txt | tr -s '\t' '\n' | sort -n | uniq | wc -l` = 150 were dropped at least once  
Of the 388 loci, `grep -c "^$" ../treeshrink/output_ts.txt` = 102 had no drops  
Use `awk '{ print NF }' ../treeshrink/output_ts.txt` to tabulate drops per line; on average there were 3.2 drops per locus with a drop (std dev = 1.8; max = 8)  
All three of the outgroups were removed from five of the loci  

After removing samples from sequences, align and run phylogenetic analysis again, this time including likelihood mapping  
Make a folder for the run `run2/` and link the `temp_input` folder there  
```s
qsub -v align_dir="temp_input",realign="y",clean="y",analysis="a",poly="y",likemap="y" /path/to/scripts/align_phylo.sh
```

The filtered alignment consisted of 388 loci with a total length of 582,030 bp and 12.4% missing data  

Plot as before  

### Likelihood mapping
As for the backbone, assess the phylogenetic signal for loci using the results from likelihood mapping  

Run from the `run2/` folder that has a `likemap` folder  
```s
ls in_align/ | cut -f 1 -d "_" > loci.txt
echo -e "Locus\tCount1\tCount2\tCount3\tCount4\tCount5\tCount6\tCount7" > likemap_counts.txt
for locus in $(cat loci.txt); do
counts=$(cut -f 8 likemap/lm_"$locus"*.lmap.quartetlh | sort -n | uniq -cd | tr -s ' ' | cut -f 2 -d ' ' | tr '\n' ' ')
paste <(echo $locus) <(echo $counts | tr ' ' '\t') >> likemap_counts.txt
done
```

Only 14 loci had > 40% of quartets in area 7, suggesting most loci were phylogenetically informative  
On average, loci had 20.2% of quartets in area 7 (std dev = 9.9; max = 63.8; min = 0.8)  
No loci had > 11% of quartets in the partly resolved areas (4, 5 and 6)  

Again, there doesn't seem to be a good reason to drop loci from the analysis under suspicion they are lowering concordance because of phylogenetic noise or lack of resolution  


# Distances
To visualise distances between samples prior to phasing, sequence alignments can be generated using all original samples and filtered loci  
The alignments can be used to construct a distance matrix for input to SplitsTree4 v. 4.17.1 and generation of a network  

Since "hybrids" are being included, use the consensus sequences from HybPhaser for calculating distances and including heterozygosity  
Note: dist.dna from the ape package treats ambiguity codes as missing data; this may be ignored (use dist.dna) or addressed by solving ambiguous positions randomly or using a different method (e.g. pofadinr)  

In the `combined_secapr` folder, make a folder `sequences_dist` to collect the sequences for alignment and cleaning  

For the A353 filtered loci, link sequences (**consensus**) keeping only the samples that are not outgroups
```s
grep -vf outgroup.txt samples.txt > grab_samples.txt
for sample in $(cat grab_samples.txt); do
ln -s /path/to/A353_secapr/HybPhaser/03_sequence_lists_filtered/samples_consensus/"$sample"_* sequences_dist/"$sample"_A353.fasta
done
```

Similarly for the OzBaits loci (using `<(grep ^3 samples.txt)` and appropriate paths), but remove the overlapping loci afterward from `sequences_dist`  
```s
python3 /path/to/scripts/remove_fastas.py -f ../shared_loci.txt *ozbaits.fasta
for file in mod_*; do
mv "$file" "${file/mod_/}"
done
```

Make a folder `distances` in the `combined_secapr` folder to run alignment and locus trees  
Make a folder `temp_input` in `distances`  
Run a script to collate the samples and loci run from the `temp_input` folder
```s
python3 /path/to/scripts/collate_loci_samples.py -s ../../samples.txt ../../sequences_dist/*{ozbaits,A353}.fasta
```

This results in 388 locus files for alignment and analysis, each with up to 135 samples  

Run a first analysis to align, clean and infer locus trees  
```s 
qsub -v align_dir="temp_input",realign="y",clean="y",analysis="l" /path/to/scripts/align_phylo.sh
```

Assess branch lengths for the loci using TreeShrink  
```s
singularity exec -H "$(pwd)" /path/to/phylo.sif run_treeshrink.py -t loci.treefile -m per-gene -q 0.10 -O output_ts -o treeshrink
```

Use the `temp_input` collated sequences (prior to alignment) but now drop TreeShrink outlier samples from the respective alignments  
```s
ls ../in_align/ | cut -f 1 -d "_" > loci.txt		# need to make sure the order is the same as in the first analysis
index=1
for locus in $(cat loci.txt); do
remove_line=$(sed -n "${index}p; $((index + 1))q" ../treeshrink/output_ts.txt | tr -s '\t' ',' | sed 's/,$//')
if [ ! -z "$remove_line" ]; then
python3 /path/to/scripts/remove_fastas.py "$remove_line" "${locus}.fasta"
else
echo "locus $locus does not have taxa to remove"
fi
((index+=1))
done
for file in mod_*; do
mv "$file" "${file/mod_/}"
done
```

Of the 135 samples, `cat ../treeshrink/output_ts.txt | tr -s '\t' '\n' | sort -n | uniq | wc -l` = 134 were dropped at least once  
Of the 388 loci, `grep -c "^$" ../treeshrink/output_ts.txt` = 118 had no drops  
Use `awk '{ print NF }' ../treeshrink/output_ts.txt` to tabulate drops per line; on average there were 3.0 drops per locus with a drop (std dev = 1.7; max = 9)  
All three of the outgroups were removed from three of the loci  

Run the alignment and start phylogenetic analysis (but cut off before completing by setting time to 3 hours)  
Note: remove the old `in_align` folder before running the job  
```s
qsub -l walltime=03:00:00 -v align_dir="temp_input",realign="y",clean="y",analysis="l" /path/to/scripts/align_phylo.sh
```

Using the locus alignments, generate distance matrices with pofadinr methods that account for ambiguity codes and also using the ape dist.dna function with the F84 model (resolving ambiguity codes randomly) for comparison; then calculate the average distance between samples for all loci  
(from folder `distances`; put sampleIDs and names to substitute into a `samples.tab` file; rename each `dist_out.nex` output before running subsequent commands or the file will be overwritten)  
```s
# run with the GENPOFAD method (uses ambiguity codes)
Rscript ~/scripts/align_to_distance.R -s samples.tab -p "g" in_align/*.fasta
# run with the MATCHSTATES method (uses ambiguity codes)
Rscript ~/scripts/align_to_distance.R -s samples.tab -p "m" in_align/*.fasta
# run for ape dist.dna and the F84 model (randomly resolve ambiguities)
Rscript ~/scripts/align_to_distance.R -s samples.tab -m "F84" -a "y" in_align/*.fasta
# run for ape dist.dna and the F84 model (ignore ambiguities)
Rscript ~/scripts/align_to_distance.R -s samples.tab -m "F84" in_align/*.fasta
```

The resulting matrices can be used to construct networks in SplitsTree4  


# "Hybrid" origins
The "hybrid" signal for samples in the datasets may result from actual recent hybridisation with another sample or close relative, or it could be caused by contamination with other samples in the dataset (at sampling, extraction, library preparation and/or sequencing stages)  

Determining which of these is the case is challenging, as well as which specific samples (rather than clades) may be "parents" or sources of contamination  

## Read depth
The hybrid signal may possibly result from low overall read depth and the presence of stochastic errors or minor contaminants that are not effectively eliminated in some samples as a result  

Assess whether mapping reads back to the assembled targets shows a pattern of greater heterozygosity in samples with lower read depth  
Create a `calls_mapping.txt` file in the SECAPR-like assemblies that specifies the references for each sample and the reads (adjust paths for A353 vs. OzBaits)
```s
script="/path/to/scripts/mapping.sh"
qcpath="/path/to/qc"		# make /path/to/qc/ozbaits for the OzBaits data
for sample in $(grep ^3 samples.txt); do
echo -e "${script} ${qcpath}/${sample}/${sample}_R1.fastq.gz ${qcpath}/${sample}/${sample}_R2.fastq.gz $(pwd)/stitch_targets/${sample}_targetcons.fasta" >> calls_mapping.txt
done
# the following doesn't apply to OzBaits
for sample in $(grep ^[7,8] samples.txt); do
echo -e "${script} ${qcpath}/stage1/${sample}/${sample}_R1.fastq.gz ${qcpath}/stage1/${sample}/${sample}_R2.fastq.gz $(pwd)/stitch_targets/${sample}_targetcons.fasta" >> calls_mapping.txt
done
```

Create a new folder (`mapping`) and launch the jobs  
```s
mkdir mapping && cd mapping
qsub -l ncpus=48,mem=96GB,jobfs=100GB,walltime=04:00:00,storage=gdata/nm31,wd -v calls_file="../calls_mapping.txt",cores_per="8",timeout="900" /path/to/scripts/launch_parallel.sh
```

Assess heterozygosity and read depth with samtools v. 1.17  
Use settings: minimum read depth of 12, reliable calls with at least 90% of bases consistent (e.g. as the main or alternative bases), and ambiguities scored when they are at least half the prevalance of the major base (e.g. 8 A, 4 T = W)  
(from each `mapping` folder; these commands can be put in a `samjob.sh` file and launched with qsub on Gadi)  
```s
#!/bin/bash
container="/path/to/samtools.sif"
module load singularity python3/3.10.0
for sample in $(cat ../samples.txt); do
cd "$sample"
singularity exec -H "$(pwd)" "$container" samtools consensus --mode simple \
	--ambig --min-depth 12 --call-fract 0.9 --het-fract 0.5 -o "$sample"_consensus.fasta map_sorted.bam
singularity exec -H "$(pwd)" "$container" samtools coverage map_sorted.bam >> "$sample"_coverage_file.txt
cd ..
done
python3 /path/to/scripts/count_ambig.py */*_consensus.fasta > ambig_file.txt
```
(if launching on Gadi: `qsub -q normal -l ncpus=4,mem=16GB,walltime=04:00:00,storage="gdata/nm31",wd samjob.sh`)  

To collate the depth results for comparison in a spreadsheet:
```s
sample=$(head -n 1 ../samples.txt)
head -n 1 "$sample"/"$sample"_coverage_file.txt > coverage_file.txt
for sample in $(cat ../samples.txt); do
tail -n +2 "$sample"/"$sample"_coverage_file.txt >> coverage_file.txt 
done
```

There is no apparent correlation between percent ambiguous bases and read depth, suggesting the hybrids are not the result of poor recovery and read depth uncertainty  
Contamination followed by high read depth of both sources of DNA is still a possibility  

## High copy targets
One possible way to assess "hybrids" is to attempt to assemble off-target sequences of the plastome and/or the ribosomal cistron  

Plastome sequences should not show heterozygosity in real hybrids (assuming typical biological inheritance and uniformity in cells), but if the sample is contaminated, then divergent plastomes should lead to detectable heterozygosity  
This may not be detectable if there is very little divergence between contaminants (e.g. if plastome variation in the group is small)  

The nuclear ribosomal cistron includes conserved and variable regions that may allow differentiation between putative parents  
Cistron sequences may undergo concerted evolution, so caution is needed for interpreting the presence/absence of variation  
Identical cistrons with copies in hybrids may be contamination (or real hybridity), but variation between hybrids and putative parents may point to real hybridity rather than contamination from that sample (M.D. Barrett, pers. comm.)  

Use SPAdes v. 3.15.5 to assemble reads from each sample (combining bait sets where possible) to form the base assembly from which to extract the plastome and cistron sequences  

In a new `SPAdes` directory, create a `calls_spades.txt` file with a line per sample pointing to the sets of A353 and OzBaits reads (or just A353 for the 9 samples from stage 1)  
Copy over the `samples.txt` file from the A353 dataset (all samples)  
```s
cd SPAdes
script="/path/to/scripts/spades.sh"
qc_dir="/path/to/qc"
for sample in $(grep ^[7,8] samples.txt); do
echo -e "${script} ${qc_dir}/stage1/${sample}/${sample}_R1.fastq.gz;${qc_dir}/stage1/${sample}/${sample}_R2.fastq.gz" >> calls_spades.txt
done
for sample in $(grep ^3 samples.txt); do
echo -e "${script} ${qc_dir}/${sample}/${sample}_R1.fastq.gz;${qc_dir}/${sample}/${sample}_R2.fastq.gz;${qc_dir}/ozbaits/${sample}/${sample}_R1.fastq.gz;${qc_dir}/ozbaits/${sample}/${sample}_R2.fastq.gz" >> calls_spades.txt
done
```

Launch the jobs using a script to parallelise the tasks  
```s
qsub -l ncpus=192,mem=768GB,jobfs=240GB,walltime=12:00:00,storage=gdata/nm31,wd -v calls_file="calls_spades.txt",cores_per="16",timeout="3000" /path/to/scripts/launch_parallel.sh
```

### Plastome
We can attempt to extract full plastomes from the assemblies using GetOrganelle v. 1.7.7.0 (in a Singularity container)  
Make a new directory `plastome`  
Copy over the `samples.txt` file from the A353 dataset (all samples)  
Note that GetOrganelle requires reference databases that need to be downloaded first (use copyq on Gadi)  
Put the following in a simple shell script (`getdb.sh`)
```s
#!/bin/bash
container="/path/to/getorganelle.sif"
module load singularity
singularity exec -H "$(pwd)" "$container" get_organelle_config.py --config-dir "$(pwd)"/.GetOrganelle \
	--add embplant_pt,embplant_mt
```
Launch with `qsub -q copyq -l ncpus=1,mem=8GB,walltime=01:00:00,storage="gdata/nm31",wd getdb.sh`  

Once the databases are downloaded, create a `calls.txt` file to run GetOrganelle for each SPAdes assembly folder
```s
script="/path/to/scripts/plastome_from_assembly.sh"
spades_path="/path/to/SPAdes"		# where the assemblies are
for sample in $(cat $samples.txt); do
echo -e "${script} ${spades_path}/${sample} ${sample}" >> calls.txt
done
```

Launch the jobs using a script to parallelise the tasks  
```s
qsub -l ncpus=96,mem=192GB,jobfs=48GB,walltime=12:00:00,storage=gdata/nm31,wd -v calls_file="calls.txt",cores_per="4" /path/to/scripts/launch_parallel.sh
```
Most assembled in less than a minute, but one sample timed-out and could not be properly recovered (376681); the resulting graph had 199 nodes that hit the reference *Theobroma cacao* plastome (NC_014676.2), totalling ~100 kb  

The files `initial_assembly_graph.fastg` in each `getorganelle_out` assembly folder can be removed (these are duplicating the SPAdes files) to save space  

Evaluate the completeness of the assemblies
```s
for sample in $(cat samples.txt); do
if [ -f $sample/getorganelle_out/*1.path_sequence.fasta ]; then
python3 ~/scripts/seq_stats.py $sample/getorganelle_out/*path_sequence.fasta | grep "Total" >> summary.txt
else
echo -e "\nNo assembly for sample $sample\n"
fi
done
```
(only 376681 lacked a recovered assembly)  

Plastome recovery was highly variable across samples, ranging from 1.4 kb to 162 kb (average 42 kb, std dev 44.5 kb)  
25 samples had > 100 kb recovered  
Only one sample had two paths recovered (orientations): 80527 (largely complete 162 kb)  

To evaluate heterozygosity, map reads for each sample against the assembly for that sample  
Note: this may underestimate heterozygosity in cases where the assembly consists of multiple copies of regions of the chloroplast as independent contigs  
Another approach might be to map to a reference instead (e.g. *Theobroma cacao*)  

Map A353 and OzBaits separately, then combined (temporarily concatenate the two sets of reads into one for the job, then delete afterward)  
```s
script="/path/to/scripts/mapping.sh"
qcpath="/path/to/qc"
mkdir temp_reads
for sample in $(grep ^3 samples.txt); do
echo -e "$script ${qcpath}/${sample}/${sample}_R1.fastq.gz ${qcpath}/${sample}/${sample}_R2.fastq.gz $(pwd)/${sample}/getorganelle_out/*1.path_sequence.fasta" >> calls_mapping_A353.txt
echo -e "$script ${qcpath}/ozbaits/${sample}/${sample}_R1.fastq.gz ${qcpath}/ozbaits/${sample}/${sample}_R2.fastq.gz $(pwd)/${sample}/getorganelle_out/*1.path_sequence.fasta" >> calls_mapping_ozbaits.txt
cat ${qcpath}/${sample}/${sample}_R1.fastq.gz ${qcpath}/ozbaits/${sample}/${sample}_R1.fastq.gz > temp_reads/${sample}_R1.fastq.gz
cat ${qcpath}/${sample}/${sample}_R2.fastq.gz ${qcpath}/ozbaits/${sample}/${sample}_R2.fastq.gz > temp_reads/${sample}_R2.fastq.gz
echo -e "$script $(pwd)/temp_reads/${sample}_R1.fastq.gz $(pwd)/temp_reads/${sample}_R2.fastq.gz $(pwd)/${sample}/getorganelle_out/*1.path_sequence.fasta" >> calls_mapping_combined.txt
done
for sample in $(grep ^[7,8] samples.txt); do
echo -e "$script ${qcpath}/stage1/${sample}/${sample}_R1.fastq.gz ${qcpath}/stage1/${sample}/${sample}_R2.fastq.gz $(pwd)/${sample}/getorganelle_out/*1.path_sequence.fasta" >> calls_mapping_A353.txt
done
```
Manually remove sample 376681 (failed)  

Make directories and launch the jobs using the script to parallelise tasks  
```s
mkdir mapping_A353
mkdir mapping_ozbaits
mkdir mapping_combined
cd mapping_A353 && qsub -l ncpus=48,mem=96GB,jobfs=100GB,walltime=04:00:00,storage=gdata/nm31,wd -v calls_file="../calls_mapping_A353.txt",cores_per="8",timeout="900" /path/to/scripts/launch_parallel.sh && cd ..
cd mapping_ozbaits && qsub -l ncpus=48,mem=96GB,jobfs=100GB,walltime=04:00:00,storage=gdata/nm31,wd -v calls_file="../calls_mapping_ozbaits.txt",cores_per="8",timeout="900" /path/to/scripts/launch_parallel.sh && cd ..
cd mapping_combined && qsub -l ncpus=48,mem=96GB,jobfs=100GB,walltime=04:00:00,storage=gdata/nm31,wd -v calls_file="../calls_mapping_combined.txt",cores_per="8",timeout="900" /path/to/scripts/launch_parallel.sh && cd ..
```

Assess read depth and heterozygosity  
(from each `mapping` folder; these commands can be put in a `samjob.sh` file and launched with qsub on Gadi)  
```s
#!/bin/bash
container="/path/to/samtools.sif"
module load singularity python3/3.10.0
for sample in $(cat ../samples.txt); do
	if [ -d "$sample" ]; then
		cd "$sample"
		singularity exec -H "$(pwd)" "$container" samtools consensus --mode simple \
			--ambig --min-depth 12 --call-fract 0.9 --het-fract 0.5 -o "$sample"_consensus.fasta map_sorted.bam
		singularity exec -H "$(pwd)" "$container" samtools coverage map_sorted.bam >> "$sample"_coverage_file.txt
		cd ..
	fi
done
python3 /path/to/scripts/count_ambig.py */*_consensus.fasta > ambig_file.txt
```
(if launching on Gadi: `qsub -q normal -l ncpus=4,mem=16GB,walltime=04:00:00,storage="gdata/nm31",wd samjob.sh`)  

To collate the depth results for comparison in a spreadsheet:
```s
sample=$(grep ^3 ../samples.txt | head -n 1)
paste <(echo "sample") <(head -n 1 "$sample"/"$sample"_coverage_file.txt) > coverage_file.txt
for sample in $(cat ../samples.txt); do		# use $(grep ^3 ../samples.txt) for OzBaits and combined
tail -n +2 "$sample"/"$sample"_coverage_file.txt | while read line; do echo -e "$sample\t$line" >> coverage_file.txt; done 
done
```

Overall, samples with some of the highest percentage ambiguities (including T. stelligera) only recovered a few kb of plastome sequence, which may suggest that the recovered contigs are not legitimate plastome (e.g. captures or off-target regions)  
The three highest percentage ambiguities are "hybrid" samples, but all had plastome recovery < 3 kb  

There is a trend toward fewer ambiguities as length of recovery increases, consistent with ambiguities at low recovery potentially reflecting illegitimate plastome sequences (e.g. poor enrichment) rather than true plastome heterozygosity  

Given the variability in recovery and the uncertainty around detecting outliers, it is hard to conclude that the heterozygosity clearly indicates contamination for most samples, though some are suspicious (e.g. T. stelligera)  

As a side note, there were apparent differences in plastome coverage for the two bait sets, with some samples doing better for one or the other  

### Ribosomal cistron
Use a *Theobroma cacao* cistron sequence (GenBank accession JQ228369.1) as a reference for extracting the nuclear ribosomal sequence contigs from the SPAdes assemblies  

Use a portion of the sequence starting from the start of 18S (coordinate 158) and extending to the end of 26S (coordinate 5924) to create a fasta reference  
```s
python3 ~/scripts/genbank_parse.py -c 158..5924 Theobroma_JQ228369.gb
```
Manually rename the file `Theo_ribo.fasta` and change the fasta header to `>Theobroma_cacao`  
The resulting reference is 5767 bp long  

Create a folder `ribosomal` and move the reference there  
Copy over the `samples.txt` file from the A353 dataset (all samples)  

In a folder `ribosomal/blast`, blast the SPAdes assemblies with the reference and record contigs that hit, keeping hits > 90% identity and longer than 500 bp  
The following can be put in a `blast_job.sh` for submission on Gadi  
```s
#!/bin/bash
module load blast/2.11.0
for sample in $(cat ../samples.txt); do
	makeblastdb -dbtype nucl -in /path/to/SPAdes/${sample}/scaffolds.fasta -out dbnucl > /dev/null
	echo "sacc slen pident length qstart qend sstart send" | tr ' ' '\t' > ${sample}_blast_out.tab
	blastn -query ../Theo_ribo.fasta -db dbnucl -outfmt "6 sacc slen pident length qstart qend sstart send" \
	-max_target_seqs 10 | awk '$3>90' | awk '$4>500' >> ${sample}_blast_out.tab
done
rm dbnucl*
```
Launch it on Gadi: `qsub -q normal -l ncpus=1,mem=4GB,walltime=02:00:00,storage="gdata/nm31",wd blast_job.sh`

In some cases, a single contig will be hit  
In others, there may be multiple contigs hitting (broken assembly)  

Create lists of samples for these two cases:
```s
for sample in $(cat samples.txt); do
hits=$(tail -n +2 blast/${sample}_blast_out.tab | cut -f 1 | sort | uniq | wc -l)
if [ $hits -gt 1 ]; then
echo $sample >> samples_stitch.txt
else
echo $sample >> samples_single.txt
fi
done
```
There are 88 samples with single hit contigs and 65 samples with multiple  

Use the 88 samples for a more reliable starting set of sequences to align  

Use the blast results to extract the single contig for each sample into a new `extracts` folder  
Note: the stitching script will only pull out a single contig if there is only one hit  
```s
mkdir extracts && cd extracts
for sample in $(cat ../samples_single.txt); do
ln -s ../blast/${sample}_blast_out.tab .
ln -s /path/to/SPAdes/${sample}/scaffolds.fasta ${sample}.fa
done
python3 /path/to/scripts/stitch_contigs.py -s ../samples_single.txt *.fa
```
The links can be removed afterward (`rm *.fa *.tab`)  

To ensure contigs are oriented to match the reference and start near the same location, run a script using the *Theobroma* reference (`Theo_ribo.fasta`) to adjust them (starting ~200 bp before the start of 18S)  
To make alignment less ambiguous, also limit the contigs to not extend much beyond the end of 26S (~100 bp)  
To do this, create sub-references for 600 bp portions of the reference at the start and end (plus a 26S start piece in case of poor samples)  
```s
python3 ~/scripts/genbank_parse.py -c 158..757 Theobroma_JQ228369.gb && mv Theobroma_extract.fasta Theo_ribo_18S_start.fasta
python3 ~/scripts/genbank_parse.py -c 5325..5924 Theobroma_JQ228369.gb && mv Theobroma_extract.fasta Theo_ribo_26S_end.fasta
python3 ~/scripts/genbank_parse.py -c 2673..3272 Theobroma_JQ228369.gb && mv Theobroma_extract.fasta Theo_ribo_26S_start.fasta
```

From a new folder `ribosomal/trimmed`, put the following in a bash script `trim_job.sh`  
This relies on the Python scripts `reform_contig.py` and `fasta_extract.py`  
```s
#!/bin/bash
module load python3/3.10.0
for sample in $(cat ../samples_single.txt); do
	contig="../extracts/${sample}_stitched.fasta"
	type="normal"
	makeblastdb -dbtype nucl -in $contig -out dbnucl > /dev/null
	blastn -query ../Theo_ribo_18S_start.fasta -db dbnucl -outfmt "6 sstart send" | head -n 1 > temp
	if [ ! -s temp ]; then		# if there was no hit to the start of 18S
		echo "Sample $sample has no hit to the start of 18S! Using 26S..."
		blastn -query ../Theo_ribo_26S_start.fasta -db dbnucl -outfmt "6 sstart send" | head -n 1 > temp
		type="poor"
	fi
	if [ ! -s temp ]; then		# if there was still no hit
		echo "Sample $sample has no hit to the start of 26S either! Moving on..."
		continue
	fi
	sstart=$(cut -f 1 temp)
	send=$(cut -f 2 temp)
	if [ $sstart -gt $send ]; then		# reverse orientation (shouldn't happen)
		python3 /path/to/scripts/reform_contig.py $contig 1 yes > /dev/null
		makeblastdb -dbtype nucl -in new_contig.fasta -out dbnucl > /dev/null
		if [ $type == "normal" ]; then
			blastn -query ../Theo_ribo_18S_start.fasta -db dbnucl -outfmt "6 sstart send" | head -n 1 > temp
		else
			blastn -query ../Theo_ribo_26S_start.fasta -db dbnucl -outfmt "6 sstart send" | head -n 1 > temp
		fi
		sstart=$(cut -f 1 temp)
		contig="new_contig.fasta"
	fi
	if [ $sstart -gt 200 ] && [ $type == "normal" ]; then	# need to reset start
		python3 /path/to/scripts/reform_contig.py $contig $((sstart - 200)) > /dev/null
		makeblastdb -dbtype nucl -in new_contig.fasta -out dbnucl > /dev/null
		contig="new_contig.fasta"
	fi
	len=$(tail -n +2 $contig | tr '\n' ' ' | sed 's/[[:space:]]//g' | wc -m)
	blastn -query ../Theo_ribo_26S_end.fasta -db dbnucl -outfmt "6 sstart send" | head -n 1 > temp2
	if [ ! -s temp2 ]; then		# if there was no hit to the end of 26S
		end=$len
	else
		send=$(cut -f 2 temp2)
		if [ $len -gt $((send + 100)) ]; then	# contig extends further
			end=$((send + 100))
		else
			end=$len
		fi
	fi
	# now trim
	if [ $type == "normal" ]; then
		start=1
	elif [ $sstart -gt 1000 ]; then		# hit 26S more than 1 kbp after start of contig
		start=$((sstart - 1000))
	else
		start=1
	fi
	python3 /path/to/scripts/fasta_extract.py $contig -c $start..$end
	mv extract.fasta ${sample}_trimmed.fasta
	rm dbnucl* temp*
	if [ -f new_contig.fasta ]; then rm new_contig.fasta; fi
done
```
(if launching on Gadi: `qsub -q normal -l ncpus=1,mem=8GB,walltime=01:00:00,storage="gdata/nm31",wd trim_job.sh`)  

Change the fasta headers to match sample number  
```s
for file in *_trimmed.fasta; do
sed -i 's/^>sequence_.*from />/' $file
done
```

Only one sample was "poor" (didn't hit the start of 18S): 376692 L. fitzgibbonii KS 1651 (contig was > 6 kb though)  
After trimming, the 88 samples were on average 6.03 kb (std dev = 0.39 kb; min = 4.36 kb; max = 6.56 kb)  

Run a phylogenetic analysis of the aligned cistrons, but only for samples from the ingroup (exclude any *Commersonia*, *Androcalva* and *Seringia*) to improve alignment of more variable regions  

Make a `phylo` folder and concatenate the cistrons into a `temp_input` folder in it  
Use a `samples_ingroup.txt` file to select the files of interest (73 samples)  
Launch the job but only run the concatenation analysis  
```s
for sample in $(cat samples_ingroup.txt); do cat ../trimmed/${sample}_trimmed.fasta >> temp_input/ribosomal.fasta; done
qsub -l ncpus=8 -v align_dir="temp_input",realign="y",clean="y",analysis="c" /path/to/scripts/align_phylo.sh
```

The filtered alignment was 6128 bp with 2.3% missing data  
Samples had on average 2.30% missing data (std dev = 6.48%; min = 0%; max = 28.93%)  
Make an `outgroup.txt` file with the three *Hannafordia* sample IDs for plotting, and a `samples.tab` file with the IDs and names  

For an additional plot, collapse nodes with < 80% Ultrafast bootstrap support  
```s
nw_ed concat.treefile "i & b < 80" o > concat_collapsed.treefile
```

The backbone of the resulting tree was partly unresolved, but some of the same main groups are evident (though the sampling is limited)  

Make a simple consensus sequence of the alignment in the `in_align` folder to use for mapping for each sample  
```s
python3 ~/scripts/consensus.py ribosomal_clean.fasta
```

For assessing ambiguities and read depth, run mapping for all samples to the consensus of the full alignment  

Map A353 and OzBaits separately, then combined (temporarily concatenate the two sets of reads into one for the job, then delete afterward)  
```s
script="/path/to/scripts/mapping.sh"
qcpath="/path/to/qc"
mkdir temp_reads
for sample in $(grep ^3 samples.txt); do
echo -e "$script ${qcpath}/${sample}/${sample}_R1.fastq.gz ${qcpath}/${sample}/${sample}_R2.fastq.gz $(pwd)/phylo/in_align/consensus.fasta" >> calls_mapping_A353.txt
echo -e "$script ${qcpath}/ozbaits/${sample}/${sample}_R1.fastq.gz ${qcpath}/ozbaits/${sample}/${sample}_R2.fastq.gz $(pwd)/phylo/in_align/consensus.fasta" >> calls_mapping_ozbaits.txt
cat ${qcpath}/${sample}/${sample}_R1.fastq.gz ${qcpath}/ozbaits/${sample}/${sample}_R1.fastq.gz > temp_reads/${sample}_R1.fastq.gz
cat ${qcpath}/${sample}/${sample}_R2.fastq.gz ${qcpath}/ozbaits/${sample}/${sample}_R2.fastq.gz > temp_reads/${sample}_R2.fastq.gz
echo -e "$script $(pwd)/temp_reads/${sample}_R1.fastq.gz $(pwd)/temp_reads/${sample}_R2.fastq.gz $(pwd)/phylo/in_align/consensus.fasta" >> calls_mapping_combined.txt
done
for sample in $(grep ^[7,8] samples.txt); do
echo -e "$script ${qcpath}/stage1/${sample}/${sample}_R1.fastq.gz ${qcpath}/stage1/${sample}/${sample}_R2.fastq.gz $(pwd)/phylo/in_align/consensus.fasta" >> calls_mapping_A353.txt
done
```

Make directories and launch the jobs using the script to parallelise tasks  
```s
mkdir mapping_A353
mkdir mapping_ozbaits
mkdir mapping_combined
cd mapping_A353 && qsub -l ncpus=48,mem=96GB,jobfs=100GB,walltime=04:00:00,storage=gdata/nm31,wd -v calls_file="../calls_mapping_A353.txt",cores_per="8",timeout="900" /path/to/scripts/launch_parallel.sh && cd ..
cd mapping_ozbaits && qsub -l ncpus=48,mem=96GB,jobfs=100GB,walltime=04:00:00,storage=gdata/nm31,wd -v calls_file="../calls_mapping_ozbaits.txt",cores_per="8",timeout="900" /path/to/scripts/launch_parallel.sh && cd ..
cd mapping_combined && qsub -l ncpus=48,mem=96GB,jobfs=100GB,walltime=04:00:00,storage=gdata/nm31,wd -v calls_file="../calls_mapping_combined.txt",cores_per="8",timeout="900" /path/to/scripts/launch_parallel.sh && cd ..
```

Assess read depth and heterozygosity, collating results  
(from each `mapping` folder; these commands can be put in a `samjob.sh` file and launched with qsub on Gadi)  
To attempt to better represent additional copies present, lower the het fraction to report if there are base calls at least 15% of the top call  
```s
#!/bin/bash
container="/path/to/samtools.sif"
module load singularity python3/3.10.0
for sample in $(cat ../samples.txt); do
	if [ -d "$sample" ]; then
		cd "$sample"
		singularity exec -H "$(pwd)" "$container" samtools consensus --mode simple \
			--ambig --min-depth 12 --het-fract 0.15 -o "$sample"_consensus.fasta map_sorted.bam
		singularity exec -H "$(pwd)" "$container" samtools coverage map_sorted.bam >> "$sample"_coverage_file.txt
		cd ..
	fi
done
python3 /path/to/scripts/count_ambig.py */*_consensus.fasta > ambig_file.txt
sample=$(grep ^3 ../samples.txt | head -n 1)
paste <(echo "sample") <(head -n 1 "$sample"/"$sample"_coverage_file.txt) > coverage_file.txt
for sample in $(cat ../samples.txt); do
	if [ -d "$sample" ]; then
		tail -n +2 "$sample"/"$sample"_coverage_file.txt | while read line; do echo -e "$sample\t$line" >> coverage_file.txt; done
	fi
done
```
(if launching on Gadi: `qsub -q normal -l ncpus=4,mem=16GB,walltime=04:00:00,storage="gdata/nm31",wd samjob.sh`)  

There was a slight trend toward fewer ambiguities at higher read depth, with the OzBaits data having substantially higher mean read depth; this may be deceptive, as much of the locus could have low depth with a single section having massive depth  
When compared with effective mapping length (no gaps or Ns; most samples had > 5 kbp mapped), some samples had higher percentages of ambiguities (> 0.5%), but there wasn't a clear break  

For the OzBaits reads, there were 10 samples, of which four were suspected "hybrids"  
For the A353 reads, there were 26 samples, of which eight were suspected "hybrids"  

Use the mapping **consensus** files for each sample to run a larger alignment, which can be put into a phylogenetic analysis, though the fact it is consensus and there are ambiguities may mean poor resolution  
The alignment can be used to evaluate putative parental ribosomal copies  

Make a `run2` folder in the `phylo` folder and concatenate the consensus files from the different mappings into a `temp_input` folder in it  
Note: the consensus files have fasta headers "consensus" that should be changed to sampleIDs  
Make a new `samples_ingroup.txt` file in the `run2` folder, again excluding any *Androcalva*, *Commersonia* and *Seringia*  
Note: one stage1 sample had > 70% missing data and low read depth and is unreliable (79869 G. ledifolia CW 2430), so exclude it  
Launch the job but only run the concatenation analysis  
```s
for sample in $(grep ^3 samples_ingroup.txt); do cat ../../mapping_ozbaits/${sample}/${sample}_consensus.fasta | sed "s/>consensus/>${sample}_ozbaits/" >> temp_input/ribosomal.fasta; done
for sample in $(cat samples_ingroup.txt); do cat ../../mapping_A353/${sample}/${sample}_consensus.fasta | sed "s/>consensus/>${sample}_A353/" >> temp_input/ribosomal.fasta; done
qsub -l ncpus=8 -v align_dir="temp_input",realign="y",clean="y",analysis="c" /path/to/scripts/align_phylo.sh
```

The filtered alignment had 263 "taxa" and was 6,016 bp with 6.0% missing data  
A353 samples had on average 5.8% missing data (std dev = 7.0%; min = 0%; max = 43%)  
OzBaits samples had on average 6.3% missing data (std dev = 9.6%; min = 0%; max = 42%)  

Collapse branches with less than 80% Ultrafast bootstrap support  
```s
nw_ed concat.treefile "i & b < 80" o > concat_collapsed.treefile
```

### ITS
Given there is indication of multiple copies in some samples, it may be insightful to extract only the ITS region (as variable enough and short enough for overlapping read coverage and phasing) and attempt to phase copies using mapped reads  
In addition, variation in ITS may include gaps in the alignment for small numbers of samples (informative), which would be removed with alignment cleaning (as in the analysis for the full cistron above)  

Use the reference to grab 50 bp portions to delimit the edges of the ITS region (plus 5.8S)  
```s
python3 ~/scripts/genbank_parse.py -c 1908..1957 Theobroma_JQ228369.gb && mv Theobroma_extract.fasta Theo_ribo_ITS_start.fasta
python3 ~/scripts/genbank_parse.py -c 2623..2672 Theobroma_JQ228369.gb && mv Theobroma_extract.fasta Theo_ribo_ITS_end.fasta
```

Starting with the trustworthy single contig ribosomal cistron assemblies (trimmed and re-oriented; 73 ingroup samples), extract the ITS region  
A preliminary alignment suggested four samples had long insertions in ITS (possible misassemblies), so don't include them when extracting (376679, 376716, 376751, 376757)  
One sample was missing most of ITS, so also exclude it (376682)  
This leaves 68 ingroup samples (put IDs in a `ITS_ref_samples.txt` file)  

Make a new directory (`ITS`) and run blast to extract the ITS from each sample  
Note: use `-task blastn` to accommodate the short queries  
```s
#!/bin/bash
for sample in $(cat ../ITS_ref_samples.txt); do
	contig="../trimmed/${sample}_trimmed.fasta"
	makeblastdb -dbtype nucl -in $contig -out dbnucl > /dev/null
	blastn -task blastn -query ../Theo_ribo_ITS_start.fasta -db dbnucl -outfmt "6 sstart send" | head -n 1 > temp
	itsstart=$(cut -f 1 temp)
	blastn -task blastn -query ../Theo_ribo_ITS_end.fasta -db dbnucl -outfmt "6 sstart send" | head -n 1 > temp2
	itsend=$(cut -f 2 temp2)
	python3 ~/scripts/fasta_extract.py $contig -c $itsstart..$itsend
	mv extract.fasta ${sample}_ITS.fasta
	rm dbnucl* temp*
	sed -i 's/^>sequence_.*from />/' ${sample}_ITS.fasta
done
```

Align with MAFFT and generate a consensus  
```s
mafft --auto <(cat *.fasta) > temp_align.fasta
python3 ~/scripts/consensus.py temp_align.fasta
rm temp_align.fasta && mv consensus.fasta consensus_ITS.fasta
```
The resulting sequence is 851 bp, with 30 `N`s and 29 ambiguities (would be mapped as `N` with BBMap)  

Run mapping of reads to the consensus, then assemble the resulting mapped reads with SPAdes to attempt to extract ITS copies  
Use reads from both bait sets  

In a new `ITS/map_assemble` directory, create a `calls_map_assemble.txt` file with a line per sample pointing to the sets of A353 and OzBaits reads (or just A353 for the 5 samples from stage 1)  
Copy over the `samples_ingroup.txt` file (134 samples) from the `phylo/run2` folder  
(adjust identity value for BBMap to 0.8 for retrieving reads divergent from the consensus)
```s
script="/path/to/scripts/map_assemble.sh"
qc_dir="/path/to/qc"
ref="/path/to/ITS/consensus_ITS.fasta"
for sample in $(grep ^[7,8] samples_ingroup.txt); do
echo -e "${script} ${qc_dir}/stage1/${sample}/${sample}_R1.fastq.gz;${qc_dir}/stage1/${sample}/${sample}_R2.fastq.gz ${ref} 0.8" >> calls_map_assemble.txt
done
for sample in $(grep ^3 samples_ingroup.txt); do
echo -e "${script} ${qc_dir}/${sample}/${sample}_R1.fastq.gz;${qc_dir}/${sample}/${sample}_R2.fastq.gz;${qc_dir}/ozbaits/${sample}/${sample}_R1.fastq.gz;${qc_dir}/ozbaits/${sample}/${sample}_R2.fastq.gz ${ref} 0.8" >> calls_map_assemble.txt
done
```

Launch the jobs using a script to parallelise the tasks  
```s
qsub -l ncpus=48,mem=96GB,jobfs=100GB,walltime=12:00:00,storage=gdata/nm31,wd -v calls_file="calls_map_assemble.txt",cores_per="8",timeout="1000" /path/to/scripts/launch_parallel.sh
```

For 3 samples (376706, 376747, 376776), SPAdes failed to produce contigs, though there were reads mapping  
(one of those samples was the "hybrid" 376776 T. x formosa, which showed heterozygosity)  

83 samples had a single contig, 38 had 2, 9 had 3, and 1 had 4  

Extract the ITS portion in the correct orientation for alignment  

For samples with single contigs, use the scaffolds files directly  

(from the `map_assemble` folder)  
```s
for sample in $(cat samples_ingroup.txt); do
if [ $(grep ">" "$sample"/scaffolds.fasta | wc -l) -eq 1 ]; then
cp "$sample"/scaffolds.fasta temp_"$sample"_contig.fasta
sed -i "s/^>.*/>$sample/" temp_"$sample"_contig.fasta
fi
done
```

For samples with multiple contigs, check blast hits on the reference: if there is no overlap, stitch contigs together separated by 10 "N"s  
If there is complete overlap, keep the first (longest) contig  

Make a little script (`stitch.py`) to stitch contigs (arg1) together, separated by 10 "N"s in a specified order (arg2) and orientation (arg3)  
e.g. `python3 stitch.py scaffolds.fasta 2,1 plus,minus` for stitching the second contig in forward orientation before the first in the reverse  
```python
#!/usr/bin/env python3
import sys
from Bio import SeqIO
from Bio.Seq import Seq
order = sys.argv[2].split(',')
orient = sys.argv[3].split(',')
contigs = []
with open(sys.argv[1], 'r') as contig_file:
	for entry in SeqIO.parse(contig_file, 'fasta'):
		contigs.append(entry)

seqs = []
with open('temp_stitch.fasta', 'w') as outfile:
	for index, item in enumerate(order):
		if orient[index] == 'plus':
			seqs.append(contigs[int(item) - 1].seq)
		else:
			seqs.append(contigs[int(item) - 1].seq.reverse_complement())

	contigs[0].seq = Seq(10 * 'N').join(seqs)
	SeqIO.write(contigs[0], outfile, 'fasta')
```

Run the stitching for samples with 2 contigs  
Put the following in a short bash script and run as `bash two_stitch.sh`
```s
#!/bin/bash
for sample in $(cat samples_ingroup.txt); do
	if [ $(grep ">" "$sample"/scaffolds.fasta | wc -l) -eq 2 ]; then
		makeblastdb -dbtype nucl -in "$sample"/scaffolds.fasta -out dbnucl > /dev/null
		blastn -query ../consensus_ITS.fasta -db dbnucl -outfmt "6 sacc qstart qend sstart send" > temp
		hit1start=$(head -n 1 temp | cut -f 2)
		hit1end=$(head -n 1 temp | cut -f 3)
		q1start=$(head -n 1 temp | cut -f 4)
		q1end=$(head -n 1 temp | cut -f 5)
		if [ $q1start -gt $q1end ]; then
			orient1="minus"
		else
			orient1="plus"
		fi
		hit2start=$(tail -n +2 temp | cut -f 2)
		if [ -s $hit2start ]; then 	# no hit
			echo "Issue with sample $sample"
			continue
		fi
		hit2end=$(tail -n +2 temp | cut -f 3)
		q2start=$(tail -n +2 temp | cut -f 4)
		q2end=$(tail -n +2 temp | cut -f 5)
		if [ $q2start -gt $q2end ]; then
			orient2="minus"
		else
			orient2="plus"
		fi
		if [ $hit1end -lt $hit2start ]; then	# no overlap
			python3 stitch.py "$sample"/scaffolds.fasta 1,2 "$orient1","$orient2"
			mv temp_stitch.fasta temp_"$sample"_contig.fasta
			sed -i "s/^>.*/>$sample/" temp_"$sample"_contig.fasta
			echo "$sample stitched"
		elif [ $hit1start -gt $hit2end ]; then		# no overlap
			python3 stitch.py "$sample"/scaffolds.fasta 2,1 "$orient2","$orient1"
			mv temp_stitch.fasta temp_"$sample"_contig.fasta
			sed -i "s/^>.*/>$sample/" temp_"$sample"_contig.fasta
			echo "$sample stitched"
		elif [ $hit1start -le $hit2start ] && [ $hit1end -ge $hit2end ]; then	# complete overlap
			python3 stitch.py "$sample"/scaffolds.fasta 1 "$orient1"
			mv temp_stitch.fasta temp_"$sample"_contig.fasta
			sed -i "s/^>.*/>$sample/" temp_"$sample"_contig.fasta
			echo "$sample stitched (complete overlap; only one contig kept)"
		else
			echo "$sample NOT stitched (partial overlap)"
			echo "$sample" >> overlap_samples.txt
		fi
	fi
done
```

Only one sample had an issue (376664), where the second contig didn't hit the reference  
Keep the first contig (check blast orientation: forward) by running:
```s
sample=376664
python3 stitch.py "$sample"/scaffolds.fasta 1 plus
mv temp_stitch.fasta temp_"$sample"_contig.fasta
sed -i "s/^>.*/>$sample/" temp_"$sample"_contig.fasta
``` 

Eight samples had partial overlap of the two contigs  

Examining where the contigs hit the reference in Bandage can provide insight into how to obtain a single sequence for each sample (for the 8 partially overlapping samples and the 10 with 3+ contigs)  

For most (5+ 6-; put sampleIDs in text files), keep the first contig only  
```s
for sample in $(cat temp_plus.txt); do
python3 stitch.py "$sample"/scaffolds.fasta 1 plus
mv temp_stitch.fasta temp_"$sample"_contig.fasta
sed -i "s/^>.*/>$sample/" temp_"$sample"_contig.fasta
done
for sample in $(cat temp_minus.txt); do
python3 stitch.py "$sample"/scaffolds.fasta 1 minus
mv temp_stitch.fasta temp_"$sample"_contig.fasta
sed -i "s/^>.*/>$sample/" temp_"$sample"_contig.fasta
done
```

For one (376671), the three contigs are fragmented, with the first and third half overlapping: drop the sample  
For another (376768), the first completely overlaps the third; stitch second and first; 2,1 plus,plus  
```s
sample=376768
python3 stitch.py "$sample"/scaffolds.fasta 2,1 plus,plus
mv temp_stitch.fasta temp_"$sample"_contig.fasta
sed -i "s/^>.*/>$sample/" temp_"$sample"_contig.fasta
``` 

For the remaining five samples, there was partial overlap that would require a consensus to recover most of the length of ITS  
First, stitch the contigs in order and correct orientation  
Next, split into two contigs; align them and run a consensus  

Make a simple script (`split.py`) to split a stitched contig at the point with 10 Ns
```python
#!/usr/bin/env python3
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
with open(sys.argv[1], 'r') as stitch_file, open('temp_split.fasta', 'w') as outfile:
	entry = SeqIO.read(stitch_file, 'fasta')
	seqs = entry.seq.split(10 * 'N')
	for index, item in enumerate(seqs):
		header = 'sequence_' + str(index)
		rec = SeqRecord(item, id = header, name = header, description = header)
		SeqIO.write(rec, outfile, 'fasta')
```

376703: 2,1 plus,minus  
376717: 1,2 plus,minus  
376748: 2,1 plus,minus  
376757: 2,1 minus,plus  
376769: 1,2 plus,minus  
```s
sample=...		# e.g. 376703
python3 stitch.py "$sample"/scaffolds.fasta ...		# e.g. 2,1 plus,minus
python3 split.py temp_stitch.fasta		# split the stitched at 10 Ns
mafft --auto <(cat temp_split.fasta ../consensus_ITS.fasta) > temp_align.fasta	# align split contigs with the consensus
python3 ~/scripts/remove_fastas.py "consensus" temp_align.fasta		# remove the consensus before generating a new one
mv mod_temp_align.fasta temp_align.fasta
python3 ~/scripts/consensus.py temp_align.fasta
mv consensus.fasta temp_"$sample"_contig.fasta
sed -i "s/^>.*/>$sample/" temp_"$sample"_contig.fasta
```

Run extraction of ITS for 130 samples (put the following in a script and run with `bash get_ITS.sh`)  
This will detect incorrect orientation and trim to the start and end of ITS  
It relies on the Python script `reform_contig.py`  
```s
#!/bin/bash
for sample in $(cat samples_ingroup.txt); do
	contig=temp_"$sample"_contig.fasta
	if [ -f $contig ]; then
		makeblastdb -dbtype nucl -in $contig -out dbnucl > /dev/null
		blastn -query ../consensus_ITS.fasta -db dbnucl -outfmt "6 sstart send" | head -n 1 > temp
		if [ $(cut -f 1 temp) -gt $(cut -f 2 temp) ]; then 		# wrong orientation for largest hit
			python3 ~/scripts/reform_contig.py $contig 1 yes > /dev/null
			makeblastdb -dbtype nucl -in new_contig.fasta -out dbnucl > /dev/null
			mv new_contig.fasta temp_contig.fasta
			contig=temp_contig.fasta
		fi
		blastn -task blastn -query ../../Theo_ribo_ITS_start.fasta -db dbnucl -outfmt "6 sstart send" | head -n 1 > temp1
		if [ ! -s temp1 ]; then		# if there was no hit to the start of ITS
			itsstart="1"
		else
			itsstart=$(cut -f 1 temp1)
		fi
		len=$(tail -n +2 $contig | tr '\n' ' ' | sed 's/[[:space:]]//g' | wc -m)
		blastn -task blastn -query ../../Theo_ribo_ITS_end.fasta -db dbnucl -outfmt "6 sstart send" | head -n 1 > temp2
		if [ ! -s temp2 ]; then		# if there was no hit to the end of ITS
			itsend=$len
		else
			itsend=$(cut -f 2 temp2)
		fi
		python3 ~/scripts/fasta_extract.py $contig -c $itsstart..$itsend
		mv extract.fasta ${sample}_ITS.fasta
		sed -i 's/^>sequence_.*from />/' ${sample}_ITS.fasta
	fi
done
```

Once this is done, remove the temporary files `rm dbnucl* temp*`
The ITS extractions can be aligned and checked, e.g. `mafft --auto <(cat *ITS.fasta) > temp_align.fasta`  

Check for evidence of real ITS copy heterozygosity in the mapping bam files in IGV  
If copies can be phased, duplicate samples and manually alter bases  
If errors are evident where stitching has occurred or elsewhere, correct them  

*********************
(Manual review of 130 samples... This takes a very long time and is probably not ideal, but was necessary to recover copies and avoid excessively divergent reads putatively from other regions of the genome)  
*********************

Since it became evident that hets in ITS1 and ITS2 could not be phased together given read coverage, split the ITS regions into separate sequences, and make copies for samples that have multiple  

Grab the 5.8S from *Theobroma* (put in the `ribosomal` directory)  
NOTE: after trying this, it became evident that the *Theobroma* reference is incorrectly annotated; use *Arabidopsis* (accession X52320) to correctly identify the coordinates in *Theobroma*  
```s
python3 ~/scripts/genbank_parse.py -c 497..660 Arabidopsis_X52320.gb && mv Arabidopsis_extract.fasta Arabidopsis_ribo_5S.fasta
makeblastdb -dbtype nucl -in Theobroma_JQ228369.fasta -out dbnucl > /dev/null
blastn -task blastn -query Arabidopsis_ribo_5S.fasta -db dbnucl -outfmt "6 sstart send"
# gives: 2165	2323
python3 ~/scripts/genbank_parse.py -c 2165..2323 Theobroma_JQ228369.gb && mv Theobroma_extract.fasta Theo_ribo_5S.fasta
```

Now split all the ITS sequences into ITS1 and ITS2 copies (from the `ITS/map_assemble` folder)  
(put the following in a `break_ITS.sh` script)  
```s
#!/bin/bash
for sample in $(cat samples_ingroup.txt); do
	contig="$sample"_ITS.fasta
	if [ -f $contig ]; then
		makeblastdb -dbtype nucl -in $contig -out dbnucl > /dev/null
		blastn -task blastn -query ../../Theo_ribo_5S.fasta -db dbnucl -outfmt "6 sstart send" | head -n 1 > temp
		if [ ! -s temp ]; then		# if there was no hit
			echo "Sample $sample didn't have a hit!!!!"
		else
			start=$(cut -f 1 temp)
			end=$(cut -f 2 temp)
			python3 ~/scripts/fasta_extract.py $contig -c 1..$start
			mv extract.fasta ${sample}_ITS1.fasta
			sed -i 's/^>sequence_.*from />/' ${sample}_ITS1.fasta
			len=$(tail -n +2 $contig | tr '\n' ' ' | sed 's/[[:space:]]//g' | wc -m)
			python3 ~/scripts/fasta_extract.py $contig -c $end..$len
			mv extract.fasta ${sample}_ITS2.fasta
			sed -i 's/^>sequence_.*from />/' ${sample}_ITS2.fasta
		fi
	else
		echo "Sample $sample doesn't have a contig"
	fi
done
```
The expected 4 samples didn't have contigs (376671, 376706, 376747, 376776)  
Once this is done, remove the temporary files `rm dbnucl* temp*`  

Now manually(!) change the ITS1 and ITS2 files where necessary and create copies  
For those that need two (or more) copies, change the names of the samples to add "_c1" and "_c2"   
Duplicate files and name them, e.g. ..._ITS1a.fasta ..._ITS1b.fasta  
Do the manual editing in a separate folder (`edited`)  

Align the two ITS parts separately with MAFFT  
```s
mafft --auto <(cat *ITS1*.fasta) > align_ITS1.fasta
mafft --auto <(cat *ITS2*.fasta) > align_ITS2.fasta
```

Make a directory `phylo` in the top `ITS` directory and copy the two alignments there  
```s
cp map_assemble/edited/align*.fasta phylo/
```

Run the analyses with IQ-TREE v. 2.2.2 for each alignment  
```s
iqtree -T 1 -s align_ITS1.fasta -m MFP --ufboot 1000 --runs 10 --prefix ITS1
iqtree -T 1 -s align_ITS2.fasta -m MFP --ufboot 1000 --runs 10 --prefix ITS2
```

Make an `outgroup.txt` file with the three *Hannafordia* sample IDs, and a `samples.tab` file with the IDs and names, making sure to add sample names to match the new headers in the case of multiple copies (e.g. "sampleID_c1	samplename_c1")  
Plot as before (for concatenation), adjusting file names appropriately  

For an optional additional plot, collapse nodes with < 80% Ultrafast bootstrap support  
```s
nw_ed ITS1.treefile "i & b < 80" o > ITS1_collapsed.treefile
nw_ed ITS2.treefile "i & b < 80" o > ITS2_collapsed.treefile
```

Based on the patterns in the two trees, it is possible to speculate which copies could be merged for an overall ITS analysis  
Concatenate sequences that group consistently (in clades) in the two trees  
Only keep the first copy for sequences that do not group separately  

Make a small script to join sequences (`join_ITS.py`) separated by 100 Ns, using a `merge_list.txt` file, in which the mergers (or lack thereof) are indicated in the third (ITS1) and fourth (ITS2) columns (sampleID in the first) as, e.g. 'c1,c2	c2,c1' for joining the first copy of ITS1 and the second of ITS2 and vice versa, or '-	c1' for joining the only copy of ITS1 and the first copy of ITS2, or '-	-' for joining the only copies of both  
```python
#!/usr/bin/env python3
import sys
from Bio import SeqIO
from Bio.Seq import Seq
def join_seq(sampleID, infile1, infile2, modifier):
	with open(infile1, 'r') as input1, open(infile2, 'r') as input2, open(str(sampleID) + '_combITS' + modifier + '.fasta', 'w') as outfile:
		fasta1 = SeqIO.read(input1, 'fasta')
		fasta2 = SeqIO.read(input2, 'fasta')
		fasta1.seq = fasta1.seq + Seq(100 * 'N') + fasta2.seq
		SeqIO.write(fasta1, outfile, 'fasta')

merge_dict = {
	'-': '.fasta',
	'c1': 'a.fasta',
	'c2': 'b.fasta'
}

with open(sys.argv[1], 'r') as infile:
	for line in infile:
		parts = line.strip().split()
		sampleID = parts[0]
		ITS1_merge = parts[2].split(',')
		ITS2_merge = parts[3].split(',')

		if len(ITS1_merge) > 1:
			join_seq(sampleID,
				str(sampleID) + '_ITS1' + merge_dict[ITS1_merge[0]],
				str(sampleID) + '_ITS2' + merge_dict[ITS2_merge[0]],
				'a'
				)
			if len(ITS2_merge) > 1:
				join_seq(sampleID,
					str(sampleID) + '_ITS1' + merge_dict[ITS1_merge[1]],
					str(sampleID) + '_ITS2' + merge_dict[ITS2_merge[1]],
					'b'
					)
			else:
				join_seq(sampleID,
					str(sampleID) + '_ITS1' + merge_dict[ITS1_merge[1]],
					str(sampleID) + '_ITS2' + merge_dict[ITS2_merge[0]],
					'b'
					)
		elif len(ITS2_merge) > 1:
			join_seq(sampleID,
				str(sampleID) + '_ITS1' + merge_dict[ITS1_merge[0]],
				str(sampleID) + '_ITS2' + merge_dict[ITS2_merge[0]],
				'a'
				)
			join_seq(sampleID,
				str(sampleID) + '_ITS1' + merge_dict[ITS1_merge[0]],
				str(sampleID) + '_ITS2' + merge_dict[ITS2_merge[1]],
				'b'
				)
		else:	# only one join
			join_seq(sampleID,
				str(sampleID) + '_ITS1' + merge_dict[ITS1_merge[0]],
				str(sampleID) + '_ITS2' + merge_dict[ITS2_merge[0]],
				''
				)
```

Once this is done, change the samples that have the same fasta headers for both copies to "_c1" and "_c2" (detect with: `grep ">" *combITS* | cut -f 2 -d ">" | uniq -cd`) and also make sure there aren't any single ITS files with a "_c*" in the header (detect with: `grep ">" *combITS.fasta | cut -f 2 -d ">"`)  

Now run alignment of the combined ITS  
```s
mafft --auto <(cat *combITS*.fasta) > align_combITS.fasta
```

Move the alignment to the `phylo` directory and run another analysis  
```s
cp ../map_assemble/edited/align_combITS.fasta .
iqtree -T 1 -s align_combITS.fasta -m MFP --ufboot 1000 --runs 10 --prefix combITS
```

From the combined ITS results (which are largely consistent with the separte ITS1 and ITS2 results, at least for supported groupings and evident hybrids), only some of the putative "hybrids" show clear indication of copies matching putative parents  
The ones showing hybridity are the known hybrid, L. x tepperi (parents L. baueri and L. discolor), and two others: L. sp. Desmond (parents L. compactum and T. microphylla) and T. angustifolia (parents close to T. multiflora and Lys. involucratum)  
The other hybrids either showed no/minimal ITS copy variation (L. behrii, L. compactum, L. longistamineum, T. discolor, and T. petalocalyx), failed to recover assembled contigs (T. x formosa; note: some mapping variation, but not clear enough), or showed some variation but did not clearly associate with possible parents (L. rosmarinifolium and allies)  
Note: T. discolor did not show copy variation, but the ITS assembled/recovered groups with *Lasiopetalum* samples, while the ribosomal cistron consensus sequence (see above) clustered with *Thomasia* samples  

## Loci consensus sequences
Recent hybrids should be heterozygous for each fixed difference between their parents, so comparing heterozygous positions in the hybrid consensus sequences and corresponding bases in other samples might indicate the most likely parental species or clades  

This could be done per locus  

At each heterozygous position in the "hybrid", there may be:  
(a) only one other base in other samples (singleton; no other information)  
(b) identical ambiguous bases in non-hybrid samples (shared diversity/ambiguity; still worth considering)  
(c) two groups of non-"hybrid" samples, each possessing one of the alternative bases for the ambiguity  

In case (c), the "parents" are in both groups; in case (b), there may be similarly two "parent" groups and multiple "hybrids"  

For each "hybrid", across sites within a locus and across loci, determine what samples occur the most in opposite groups (i.e. were rarely in the same "parent" group)  
Those samples are the most likely true "parents" or "sources" of the admixed sample (or closely related to them)  

First, use a combined concatenated alignment from the distance analysis (which used consensus sequences)  
Put the "hybrid" sample IDs to test in `hybrids.txt`  
Keep top combos only when they exceed 20% of all ambiguities (across sites, or across sites within loci)  
```s
# from the "distances/in_align" folder
python3 ~/scripts/combine_alignments.py -f single *.fasta
mv combine_out.fasta .. && cd ..
python3 ~/scripts/hybrid_parents.py -f hybrids.txt -t 20 combine_out.fasta > hybrid_report_concat.txt
```
For the known hybrid (L. x tepperi), the recovered most likely parental lineages match what was suspected previously, but other results are potentially not as biologically sensible  
Some "hybrids" have a more clear top parent combination, and a couple are more consistently divergent at heterozygous positions (the top count is > 50% of the ambiguities)  
The known hybrid has the greatest proportion of its ambiguities supporting the top count (71%)  

Assess individual loci, keeping top counts only when they exceed 20% of the total ambiguities, and only tally parents when there are fewer than 10 unique top parents  
```s
python3 ~/scripts/hybrid_parents.py -f hybrids.txt -t 20 -r 10 in_align/*.fasta > hybrid_report_loci.txt
```
The known hybrid shows the strongest signal as far as the percentage of loci for the top counted parent (73%), but other "hybrids" also have high top counts  
The relative difference between counts for parents is more striking in the known hybrid and a few other "hybrids", but many have relatively indistinct indications of the top parents  

Do the different bait sets produce similar results for hybrid parents?  

Assess the concatenated loci, this time only concatenating the A353 or OzBaits loci
```s
python3 ~/scripts/combine_alignments.py -f single [4-7]*.fasta
mv combine_out.fasta ../combine_out_A353.fasta
python3 ~/scripts/combine_alignments.py -f single A*.fasta
mv combine_out.fasta ../combine_out_OzBaits.fasta
cd ..
python3 ~/scripts/hybrid_parents.py -f hybrids.txt -t 20 combine_out_A353.fasta > hybrid_report_concat_A353.txt
python3 ~/scripts/hybrid_parents.py -f hybrids.txt -t 20 combine_out_OzBaits.fasta > hybrid_report_concat_OzBaits.txt
```

Assess the loci again, this time only choosing A353 or OzBaits loci
```s
python3 ~/scripts/hybrid_parents.py -f hybrids.txt -t 20 -r 10 in_align/[4-7]*.fasta > hybrid_report_loci_A353.txt
python3 ~/scripts/hybrid_parents.py -f hybrids.txt -t 20 -r 10 in_align/A*.fasta > hybrid_report_loci_OzBaits.txt
```

Though the OzBaits data have far fewer loci with > 8 heterozygous sites or passing the 20% threshold (average of 19.5 compared to A353 average of 205), they largely agree with the A353 loci in top parents  
