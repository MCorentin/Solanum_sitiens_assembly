# Solanum sitiens assembly Pipeline

Instructions and files to reproduce the de novo assembly of Solanum sitiens (accession LA1974).

## How to cite

Coming soon.

## Table of Contents

- [Solanum sitiens assembly Pipeline](#solanum-sitiens-assembly)
- [How to cite](#how-to-cite)
- [Table of Contents](#table-of-contents)
- [Data](#Data)
- [The assembly pipeline](#the-assembly-pipeline)
    - [Versions](#versions)
    - [MaSuRCA](#masurca)
    - [SSPACE](#sspace)
    - [Hybrid Scaffold](#hybrid-scaffold)
    - [pilon](#pilon)
    - [GapFiller](#gapfiller)
    - [BBMap dedupe](#bbmap-dedupe)
    - [Arcs + Link](#arcs-+-link)
- [Gene prediction with Augustus](#gene-prediction-with-augustus)
    - [Soft-masking using RepeatMasker](#soft-masking-using-repeatMasker)
    - [Running Augustus](#running-augustus)


## Data

The assembly and reads have been submitted to the Sequence Read Archive under the Bioproject:

 - "[PRJNA633104](http://www.ncbi.nlm.nih.gov/bioproject/633104 "Solanum sitiens BioProject")  Solanum sitiens cultivar:LA1974 Genome sequencing and assembly"

 - The Bionano optical mapping data is available as supplementary file of the Bioproject: "SUPPF_0000003614"

## The assembly pipeline

### Versions

- MaSuRCA             v3.2.2
- SSPACE              v1-1
- Hybrid Scaffold     v1.0
- RefAligner          v1.0
- bwa                 v0.7.15
- samtools            v1.9
- Pilon               v1.22
- GapFiller           v1.10
- BBMap               v37.72
- LongRanger          v2.2.2
- Arcs                v8.25
- Links               v1.8.5

- STAR                v2.6.0c
- RepeatMasker        open-4.0.9
- Augustus            v3.3

### MaSuRCA

Generating the initial assembly from the Illumina and PacBio reads.

````
# First generate the "assemble.sh" script from the config file
masurca -g config.txt

# Then run the assemble.sh script
sh assemble.sh
````

Output: *dedup.genome.scf.fasta*

### SSPACE

Scaffolding of the sequences with the PacBio reads.

````
perl SSPACE-LongRead.pl -c ../dedup.genome.scf.fasta -p Sequel_RSII.fasta -b ./ - t 40
````

Output: *scaffolds.fasta*

### Hybrid Scaffold

Scaffolding with Bionano optical mapping.

````
perl Solve_06082017Rel/HybridScaffold/1.0/hybridScaffold.pl -n scaffolds.fasta -b exp_refineFinal1_contigs_S.cmap -c hybridScaffold_config.xml -r /path/to/RefAligner -o ./ -f -B 1 -N 1
````

Reintegrating the unmapped scaffolds with "HybridScaffold_finish.pl" (https://github.com/i5K-KINBRE-script-share/Irys-scaffolding)

````
perl Irys-scaffolding/KSU_bioinfo_lab/map_tools/hybridScaffold_finish_fasta.pl -x ../hybrid_scaffolds/exp_refineFinal1_contigs_S_bppAdjust_cmap_scaffolds_fasta_NGScontigs_HYBRID_SCAFFOLD.xmap -s ../hybrid_scaffolds/exp_refineFinal1_contigs_S_bppAdjust_cmap_scaffolds_fasta_NGScontigs_HYBRID_SCAFFOLD.fasta -f ../../scaffolds.fasta
````

Output: *scaffolds_genome_post_HYBRID_SCAFFOLD.fasta*


### Pilon

Correcting errors and misassemblies with the Illumina reads.

````
bwa index scaffolds_genome_post_HYBRID_SCAFFOLD.fasta

bwa mem -t 40 scaffolds_genome_post_HYBRID_SCAFFOLD.fasta 1490_R1.fastq 1490_R2.fastq > 1490.sam
bwa mem -t 40 scaffolds_genome_post_HYBRID_SCAFFOLD.fasta 1494_R1.fastq 1494_R2.fastq > 1494.sam

samtools view -b 1490.sam -@ 40 > 1490.bam
samtools view -b 1494.sam -@ 40 > 1494.bam

samtools sort 1490.bam -@ 40 > 1490_sorted.bam
samtools sort 1494.bam -@ 40 > 1494_sorted.bam
````

Running pilon:
````
java -jar -Xmx750G pilon-1.22.jar --genome scaffolds_genome_post_HYBRID_SCAFFOLD.fasta --frags 1490_sorted.bam --frags 1494_sorted.bam --output corrected_sitiens_assembly.fasta --outdir ./ --changes --fix all --threads 60
````

Output: *corrected_sitiens_assembly.fasta*

### GapFiller

Filling the gaps (stretch of N's) with the Illumina paired-end reads.

````
perl ./GapFiller_v1-10_linux-x86_64/GapFiller.pl -l libraries.txt -s corrected_sitiens_assembly.fasta -m 20 -i 10 -T 20 -b sitiens_pilonRound1_gapfiller
````

Output: *sitiens_pilonRound1_gapfiller.gapfilled.final.fa*

### BBMap dedupe

Removing duplicated contigs.

````
bash dedupe.sh in=sitiens_pilonRound1_gapfiller.gapfilled.final.fa out=sitiens_pilonRound1_gapfiller.gapfilled.final.deduped.fa outd=duplicateScaffolds.fasta threads=60 storequality=f absorbrc=t touppercase=t minidentity=90 minlengthpercent=0 minoverlappercent=0 maxsubs=40000 maxedits=5000 minoverlap=1000 k=31 -eoom -Xmx300G
````

Output: *sitiens_pilonRound1_gapfiller.gapfilled.final.deduped.fa*

### Arcs + Link

Scaffolding with the Chromium 10x data.

Interleaving the chromium 10x reads with LongRanger.
*/path/to/chromium_10x/Sample_PRO1846_S1R1/* is the folder containing the paired 10x fastq files.
````
# First generate the interleaved file
./longranger/2.2.2/longranger basic --id=sitiens_10x_basic --fastqs=/path/to/chromium_10x/Sample_PRO1846_S1R1/

# Then add the barcode to the file (needed for Arcs)
gunzip -c barcoded.fastq.gz | perl -ne 'chomp;$ct++;$ct=1 if($ct>4);if($ct==1){if(/(\@\S+)\sBX\:Z\:(\S{16})/){$flag=1;$head=$1."_".$2;print "$head\n";}else{$flag=0;}}else{print"$_\n" if($flag);}' > longranger_perl.fastq
````

Aligning the interleaved chromium 10x reads to the assembly:
````
bwa mem -t 50 -p sitiens_pilonRound1_gapfiller.gapfilled.final.deduped.fa longranger_perl.fastq > sitiens_pilonRound1_gapfiller.gapfilled.final.deduped.sam

samtools view -b sitiens_pilonRound1_gapfiller.gapfilled.final.deduped.sam -@ 40 > sitiens_pilonRound1_gapfiller.gapfilled.final.deduped.bam

# Sorting the bam file by read name (necessary for Arcs):
samtools sort -n sitiens_pilonRound1_gapfiller.gapfilled.final.deduped.bam -@ 40 > /path/to/sitiens_pilonRound1_gapfiller.gapfilled.final.deduped.sortedByName.bam

echo "/path/to/sitiens_pilonRound1_gapfiller.gapfilled.final.deduped.sortedByName.bam" > alignment.fof
````

Launching Arcs:
````
arcs -f sitiens_pilonRound1_gapfiller.gapfilled.final.deduped.fa -a alignment.fof -s 98 -c 5 -l 0 -z 500 -m 50-10000 -d 0 -e 30000 -r 0.05 -v 1
````

Launching Links (*empty.fof* is an empty file).
````
# Creating the checkpoint.tsv file for Links
./arcs/Examples/makeTSVfile.py sitiens_pilonRound1_gapfiller.gapfilled.final.deduped.fa.scaff_s98_c5_l0_d0_e30000_r0.05_original.gv sitiens_pilonRound1_gapfiller.gapfilled.final.deduped.Arcs.tigpair_checkpoint.tsv sitiens_pilonRound1_gapfiller.gapfilled.final.deduped.fa

# Launching Links:
./LINKS -f sitiens_pilonRound1_gapfiller.gapfilled.final.deduped.fa -b sitiens_pilonRound1_gapfiller.gapfilled.final.deduped.Arcs -s empty.fof -k 20 -l 5 -t 2 -v 1
````

Ouput: *sitiens_pilonRound1_gapfiller.gapfilled.final.deduped.Arcs.scaffolds.fa*


## Gene prediction with Augustus

Augustus predicted genes on the scaffold assembly, with RNA-seq as evidence.

### Soft-masking using RepeatMasker

The "draft_repeats_master.v5.fasta" file was obtained from SolGenomics.
ftp://ftp.solgenomics.net/tomato_genome/repeats/draft_repeats_master.v5.fasta

````
RepeatMasker --species tomato --noisy --xsmall --lib draft_repeats_master.v5.fasta sitiens_pilonRound1_gapfiller.gapfilled.final.deduped.Arcs.scaffolds.fa
````

### Running Augustus

The hints are created with the "bam2hints.pl" script from Augustus. "sitiens_RNA_Merged.sorted.bam" comes from the
alignment of the RNA-seq to the scaffold assembly with STAR. The resulting bam files are merged with "samtools merge".
````
/path/to/Augustus/auxprogs/bam2hints/bam2hints --intronsonly --in=sitiens_RNA_Merged.sorted.bam --out=introns_2.gff
````

Runing Augustus with the hints:
````
augustus --species=tomato --extrinsicCfgFile=extrinsic.M.RM.E.W.P.tomato.cfg --hintsfile=introns_2.gff --allow_hinted_splicesites=atac --alternatives-from-evidence=on --softmasking=on sitiens_pilonRound1_gapfiller.gapfilled.final.deduped7.Arcs.scaffolds.fa.masked > augustus.hints.gff
````

Extracting the sequences from the gff:
````
/path/to/Augustus/scripts/getAnnoFasta.pl augustus.hints.gff --seqfile sitiens_pilonRound1_gapfiller.gapfilled.final.deduped7.Arcs.scaffolds.fa.masked
````

Output:
- `augustus.hints.aa`           =  Amino acid sequences
- `augustus.hints.cdsexons`     =  Coding exon positions on genome
- `augustus.hints.codingseq`    =  Coding sequences only
- `augustus.hints.mrna`         =  mRNA sequences
