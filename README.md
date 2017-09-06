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
