
# A bioinformatics pipeline based on [Ruffus](http://www.ruffus.org.uk/)

Author: Ajayi Olabode (boraton2010@gmail.com)

SNPs_Analysis is based on the [Ruffus](http://www.ruffus.org.uk/) library for writing bioinformatics pipelines. Its features include:

 * Job submission on a cluster using DRMAA (currently only tested with SLURM).
 * Job dependency calculation and checkpointing.
 * Pipeline can be displayed as a flowchart.
 * Re-running a pipeline will start from the most up-to-date stage. It will not redo previously completed tasks.

## License

3 Clause BSD License. See LICENSE.txt in source repository.

## Installation

#### External dependencies

`SNPs_Analysis` depends on the following programs and libraries:

 * [python](https://www.python.org/download/releases/2.7.5/) (version 2.7.5)
 * [DRMAA](http://www.drmaa.org/) for submitting jobs to the cluster (it uses the Python wrapper to do this). 
   You need to install your own `libdrama.so` for your local job submission system. There are versions
   available for common schedulers such as Torque/PBS, [SLURM](http://apps.man.poznan.pl/trac/slurm-drmaa) and so on.
 * [bwa](http://bio-bwa.sourceforge.net/) for aligning reads to the reference genome (version 0.7.10)
 * [gatk](https://www.broadinstitute.org/gatk/) Genome Analysis Toolkit (version 3.3-0)
 * [samtools](http://www.htslib.org/) (version 0.1.2)
 * [picard](http://broadinstitute.github.io/picard/) (version 1.127)

### Important Notes:

* The pipeline assumes no known variants are available for the Base Quality Score Recalibration step and as such bootsraps a set of SNPs and Indels to provide as input at this step. Contact me if a list of high quality known variants is available for your organism as this data can improve the quality of the output
* The current script is designed to work with Paired End data
* Due to HPC queue limits on mercer, the pipeline can run on a maximum of 25 libraries at a time

# Walk-Through Analysis pipeline Steps

|   Step 1.     | Alignment – Map to Reference |
| ------------- |:-------------:|
| Tool          | BWA MEM |
| Input         | .fastq files, reference genome |
| Output        | aligned_reads.sam*|
|               | *Intermediary file, removed from final output |
| Command       | bwa mem -M -R '@RG\tID:sample_1\tLB:sample_1\tPL:ILLUMINA\tPM:HISEQ\tSM:sample_1' ref input_1 input_2 > aligned_reads.sam|


| Step 2        | Sort SAM file by coordinate, convert to BAM |
| ------------- |:-------------:|
| Tool          | Picard Tools  |
| Input         | aligned_reads.sam |
| Output        | sorted_reads.bam* |
|               | *Intermediary file, removed from final output |
| Command       | java -jar picard.jar SortSam INPUT=aligned_reads.sam OUTPUT=sorted_reads.bam SORT_ORDER=coordinate |


| Step 3        | Collect Alignment & Insert Size Metrics |
| ------------- |:-------------:|
| Tool          | Picard Tools, R, Samtools |
| Input         | sorted_reads.bam, reference genome |
| Output        | alignment_metrics.txt, insert_metrics.txt, insert_size_histogram.pdf |
| Command       | java -jar picard.jar CollectAlignmentSummaryMetrics R=ref I=sorted_reads.bam O=alignment_metrics.txt |
|               | java -jar picard.jar CollectInsertSizeMetrics INPUT=sorted_reads.bam OUTPUT=insert_metrics.txt HISTOGRAM_FILE=insert_size_histogram.pdf|

| Step 4        | Mark MarkDuplicates |
| -------------:|:-------------:|
| Tool          | Picard Tools| 
| Input         | sorted_reads.bam |
| Output        | dedup_reads.bam, metrics.txt |
|               | *Intermediary file, removed from final output |
| Command       | java -jar picard.jar MarkDuplicates INPUT=sorted_reads.bam OUTPUT=dedup_reads.bam METRICS_FILE=metrics.txt |


| Step 5        | Build BAM Index|
| -------------:|:--------------:|
| Tool          | Picard Tools   |
| Input         | dedup_reads.bam|
| Output        | dedup_reads.bai*|
|               | *Intermediary file, removed from final output|
| Command       | java -jar picard.jar BuildBamIndex INPUT=dedup_reads.bam|


| Step 6        | Create Realignment Targets |
| -------------:|:-------------:|
| Tool          | GATK          |
| Input         | dedup_reads.bam,reference genome |
| Output        | realignment_targets.list |
|               | Notes This is the first step in a two-step process of realigning around indels |
| Command       | java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R ref -I dedup_reads.bam -o realignment_targets.list |


| Step 7        | Realign Indels |
|--------------:|:---------------|
| Tool          | GATK           |
| Input         | dedup_reads.bam, realignment_targets.list, reference genome |
| Output        | realigned_reads.bam |
|               | Notes This is the first step in a two-step process of realigning around indels |
| Command       | java -jar GenomeAnalysisTK.jar -T IndelRealigner -R ref -I dedup_reads.bam -targetIntervals realignment_targets.list -o realigned_reads.bam |


| Step 8       | Call Variants |
| ------------:|:--------------|
| Tool         | GATK          |
| Input        | realigned_reads.bam, reference genome |
| Output       | raw_variants.vcf |
|              | Notes First round of variant calling. The variants identified in this step will be filtered and provided as input for Base Quality Score Recalibration |
| Command      | java -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R ref -I realigned_reads.bam -o raw_variants.vcf |


| Step 9       | Extract SNPs & Indels |
| ------------:|:-------------:|
| Tool         | GATK          |
| Input        | raw_variants.vcf, reference genome |
| Output       | raw_indels.vcf, raw_snps.vcf |
|              | Notes This step separates SNPs and Indels so they can be processed and used independently |
| Command      | java -jar GenomeAnalysisTK.jar -T SelectVariants -R ref -V raw_variants.vcf -selectType SNP -o raw_snps.vcf |
|              | java -jar GenomeAnalysisTK.jar -T SelectVariants -R ref -V raw_variants.vcf -selectType INDEL -o raw_indels.vcf |


| Step 10      | Filter SNPs |
| ------------:|:------------|
| Tool         | GATK        |
| Input        | raw_snps.vcf, reference genome |
| Output       | filtered_snps.vcf |
|              | Notes The filtering criteria for SNPs are as follows:|
|              | QD < 2.0    |
|              | FS > 60.0   |
|              | MQ < 40.0   |
|              | MQRankSum < -12.5 |
|              | ReadPosRankSum < -8.0 |
|              | SOR > 4.0   |
|              | Note: SNPs which are ‘filtered out’ at this step will remain in the filtered_snps.vcf file, however they will be marked as ‘basic_snp_filter’, while SNPs which passed the filter will be marked as ‘PASS’ |
| Command     | java -jar GenomeAnalysisTK.jar -T VariantFiltration -R ref -V raw_snps.vcf --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filterName "basic_snp_filter" -o filtered_snps.vcf |


| Step 11     | Filter Indels |
| -----------:|:--------------|
| Tool        | GATK          |
| Input       | raw_indels.vcf, reference genome |
| Output      | filtered_indels.vcf
|             | Notes The filtering criteria for SNPs are as follows:
|             | QD < 2.0      |
|             | FS > 200.0    |
|             | ReadPosRankSum < -20.0 |
|             |SOR > 10.0     |
|             | Note: Indelss which are ‘filtered out’ at this step will remain in the filtered_indels.vcf file, however they will be marked as ‘basic_indel_filter’, while Indels which passed the filter will be marked as ‘PASS’ |
| Command     | java -jar GenomeAnalysisTK.jar -T VariantFiltration -R ref -V raw_indels.vcf --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filterName "basic_indel_filter" -o filtered_indels.vcf |


| Step 12     | Base Quality Score Recalibration (BQSR) #1 |
| -----------:|:--------------|
| Tool        | GATK          |
| Input       | realigned_reads.bam, filtered_snps.vcf, filtered_indels.vcf, reference genome |
| Output      | recal_data.table* |
|             | *Intermediary file, removed from final output |
|             | Notes: BQSR is performed twice. The second pass is optional, but is required to produce a recalibration report. |
| Command     | java -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R ref -I realigned_reads.bam -knownSites filtered_snps.vcf -knownSites filtered_indels.vcf -o recal_data.table |


| Step 13     | Base Quality Score Recalibration (BQSR) #2 |
| -----------:|:--------------|
| Tool        | GATK          |
| Input       | recal_data.table, realigned_reads.bam, filtered_snps.vcf, filtered_indels.vcf, reference genome |
| Output      | post_recal_data.table | 
|             | *Intermediary file, removed from final output |
|             | Notes The second time BQSR is run, it takes the output from the first run (recal_data.table) as input |
| Command     | java -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R ref -I realigned_reads.bam -knownSites filtered_snps.vcf -knownSites filtered_indels.vcf -BQSR recal_data.table -o post_recal_data.table |


| Step 14     | Analyze Covariates |
| -----------:|:--------------|
| Tool        | GATK          |
| Input       | recal_data.table, post_recal_data.table, reference genome |
| Output      | recalibration_plots.pdf |
|             | Notes This step produces a recalibration report based on the output from the two BQSR runs |
| Command     | java -jar GenomeAnalysisTK.jar -T AnalyzeCovariates -R ref -before recal_data.table -after post_recal_data.table -plots recalibration_plots.pdf |


| Step 15     | Apply BQSR    |
| -----------:|:--------------|
| Tool        | GATK          |
| Input       | recal_data.table, realigned_reads.bam, reference genome |
| Output      | recal_reads.bam |
|             | Notes This step applies the recalibration computed in the first BQSR step to the bam file. |
| Command     | java -jar GenomeAnalysisTK.jar -T PrintReads -R ref -I realigned_reads.bam -BQSR recal_data.table -o recal_reads.bam |


| Step 16     | Call Variants |
| -----------:|:--------------|
| Tool        | GATK          |
| Input       | recal_reads.bam, reference genome |
| Output      | raw_variants_recal.vcf|
|             |Notes Second round of variant calling performed on recalibrated bam |
| Command     | java -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R ref -I recal_reads.bam -o raw_variants_recal.vcf |


| Step 17     | Extract SNPs & Indels |
| -----------:|:--------------|
| Tool        | GATK          | 
| Input       | raw_variants_recal.vcf, reference genome |
| Output      | raw_indels_recal.vcf, raw_snps_recal.vcf |
|             | Notes This step separates SNPs and Indels so they can be processed and analyzed independently |
| Command     | java -jar GenomeAnalysisTK.jar -T SelectVariants -R ref -V raw_variants_recal.vcf -selectType SNP -o raw_snps_recal.vcf |
|             | java -jar GenomeAnalysisTK.jar -T SelectVariants -R ref -V raw_variants_recal.vcf -selectType INDEL -o raw_indels_recal.vcf |


You will need to install these dependencies yourself.

I recommend using a virtual environment:

```
cd /place/to/install
virtualenv SNPs_Analysis
source SNPs_Analysis/bin/activate
pip install -U git+https://github.com/boratonAJ/SNPs_Analysis
```

If you don't want to use a virtual environment then you can just install with pip:

```
pip install -U git+https://github.com/boratonAJ/SNPs_Analysis
```

## Worked example

The `example` directory in the source distribution contains a small dataset to illustrate the use of the pipeline.

#### Get a copy of the source distribution

```
cd /path/to/test/directory
git clone https://github.com/boratonAJ/SNPs_Analysis.git
```

#### Install `SNPs_Analysis` as described above

#### Get a reference genome.

```
cd SNPs_Analysis/example
mkdir reference
# copy your reference into this directory, or make a symbolic link
# call it reference/H37rv.fa
```

#### Tell Python where your DRMAA library is 

For example (this will depend on your local settings):

```
export DRMAA_LIBRARY_PATH=/usr/local/slurm_drmaa/1.0.7-gcc/lib/libdrmaa.so
```

#### Run `SNPs_Analysis` and ask it what it will do next

```
SNPs_Analysis -n --verbose 3
```

#### Generate a flowchart diagram

```
SNPs_Analysis --flowchart pipeline_flow.png --flowchart_format png
```

#### Run the pipeline

```
SNPs_Analysis --use_threads --log_file pipeline.log --jobs 2 --verbose 3
```

## Usage

You can get a summary of the command line arguments like so:

```
SNPs_Analysis -h
usage: SNPs_Analysis     [-h] [--verbose [VERBOSE]] [-L FILE] [-T JOBNAME]
                         [-j N] [--use_threads] [-n] [--touch_files_only]
                         [--recreate_database] [--checksum_file_name FILE]
                         [--flowchart FILE] [--key_legend_in_graph]
                         [--draw_graph_horizontally]
                         [--flowchart_format FORMAT] [--forced_tasks JOBNAME]
                         [--config CONFIG] [--jobscripts JOBSCRIPTS]
                         [--version]


optional arguments:
  -h, --help            show this help message and exit
  --config CONFIG       Pipeline configuration file in YAML format, defaults
                        to pipeline.config
  --jobscripts JOBSCRIPTS
                        Directory to store cluster job scripts created by the
                        pipeline, defaults to jobscripts
  --version             show program's version number and exit

Common options:
  --verbose [VERBOSE], -v [VERBOSE]
                        Print more verbose messages for each additional
                        verbose level.
  -L FILE, --log_file FILE
                        Name and path of log file

pipeline arguments:
  -T JOBNAME, --target_tasks JOBNAME
                        Target task(s) of pipeline.
  -j N, --jobs N        Allow N jobs (commands) to run simultaneously.
  --use_threads         Use multiple threads rather than processes. Needs
                        --jobs N with N > 1
  -n, --just_print      Don't actually run any commands; just print the
                        pipeline.
  --touch_files_only    Don't actually run any commands; just 'touch' the
                        output for each task to make them appear up to date.
  --recreate_database   Don't actually run any commands; just recreate the
                        checksum database.
  --checksum_file_name FILE
                        Path of the checksum file.
  --flowchart FILE      Don't run any commands; just print pipeline as a
                        flowchart.
  --key_legend_in_graph
                        Print out legend and key for dependency graph.
  --draw_graph_horizontally
                        Draw horizontal dependency graph.
  --flowchart_format FORMAT
                        format of dependency graph file. Can be 'svg', 'svgz',
                        'png', 'jpg', 'psd', 'tif', 'eps', 'pdf', or 'dot'.
                        Defaults to the file name extension of --flowchart
                        FILE.
  --forced_tasks JOBNAME
                        Task(s) which will be included even if they are up to
                        date.
```

## Configuration file

You must supply a configuration file for the pipeline in YAML format.

Here is an example:

```
        walltime: '10:00'
        mem: 30
        modules:
            - 'snpeff/default'

# The Human Genome in FASTA format.

ref_grch37: /usr/people/ajayi/test/complexo_pipeline/example/reference/HumanTest500k_g1k_H37Rv_decoy.fasta

index_file: /usr/people/ajayi/test/complexo_pipeline/example/reference/*.nix

# The input FASTQ files.

fastqs:
   - /cip0/research/scratch/ajayi/Input_fasta_files/H37Rv1117_R1.fastq.gz
   - /cip0/research/scratch/ajayi/Input_fasta_files/H37Rv1117_R2.fastq.gz
   - /cip0/research/scratch/ajayi/Input_fasta_files/H37Rv1118_R1.fastq.gz
   - /cip0/research/scratch/ajayi/Input_fasta_files/H37Rv1118_R2.fastq.gz
   - /cip0/research/scratch/ajayi/Input_fasta_files/H37Rv1119_R1.fastq.gz
   - /cip0/research/scratch/ajayi/Input_fasta_files/H37Rv1119_R2.fastq.gz
   - /cip0/research/scratch/ajayi/Input_fasta_files/H37Rv1120_R1.fastq.gz
   - /cip0/research/scratch/ajayi/Input_fasta_files/H37Rv1120_R2.fastq.gz
   - /cip0/research/scratch/ajayi/Input_fasta_files/H37Rv1121_R1.fastq.gz
   - /cip0/research/scratch/ajayi/Input_fasta_files/H37Rv1121_R2.fastq.gz


pipeline_id: 'H37Rv'

```
