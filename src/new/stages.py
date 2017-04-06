'''
Individual stages of the pipeline implemented as functions from
input files to output files.

The run_stage function knows everything about submitting jobs and, given
the state parameter, has full access to the state of the pipeline, such 
as config, options, DRMAA and the logger.
'''

from utils import safe_make_dir
from runner import run_stage
import os
import logging
#log = logging.getLogger( __name__ )

PICARD_JAR='$PICARD_HOME/picard.jar'
GATK_JAR = '$GATK_HOME/GenomeAnalysisTK.jar'
#TRIM_PATH ='/usr/people/ajayi/test/complexo_pipeline/example/fastqs'

def java_command(jar_path, mem_in_gb, command_args):
    '''Build a string for running a java command'''
    # Bit of room between Java's max heap memory and what was requested.
    # Allows for other Java memory usage, such as stack.
    java_mem = mem_in_gb - 1
    return 'java -jar -Xmx{mem}g {jar_path} {command_args}'.format(
        jar_path=jar_path, mem=java_mem, command_args=command_args)


def run_java(state, stage, jar_path, mem, args):
    command = java_command(jar_path, mem, args)
    run_stage(state, stage, command)


class Stages(object):
    def __init__(self, state):
        self.state = state
        self.reference = self.get_options('ref_grch37')
        self.dbsnp_grch37 = self.get_options('dbsnp_grch37')
        self.mills_grch37 = self.get_options('mills_grch37')
        self.one_k_g_snps = self.get_options('one_k_g_snps')
        self.one_k_g_indels = self.get_options('one_k_g_indels')
        self.one_k_g_highconf_snps = self.get_options('one_k_g_highconf_snps')
        self.hapmap = self.get_options('hapmap')
        self.interval_grch37 = self.get_options('interval_grch37')
        self.CEU_mergeGvcf = self.get_options('CEU_mergeGvcf')
        self.GBR_mergeGvcf = self.get_options('GBR_mergeGvcf')
        self.FIN_mergeGvcf = self.get_options('FIN_mergeGvcf')

    def run_picard(self, stage, args):
        mem = int(self.state.config.get_stage_options(stage, 'mem'))
        return run_java(self.state, stage, PICARD_JAR, mem, args)

    def run_gatk(self, stage, args):
        mem = int(self.state.config.get_stage_options(stage, 'mem'))
        return run_java(self.state, stage, GATK_JAR, mem, args)


    def get_stage_options(self, stage, *options):
        return self.state.config.get_stage_options(stage, *options)

    def get_options(self, *options):
        return self.state.config.get_options(*options)


    def original_fastqs(self, output):
        '''Original fastq files'''
        pass

    def trim_galore(self, input_files, output_paired_files):
	'''Trim the paired end fastq files to the get quality data'''
	out_R1, out_R2 = output_paired_files
	read_in_R1, read_in_R2 = input_files
    	if len(input_files) != 2:
        	raise Exception("One of read pairs %s missing" % (input_files,))
    	cores = self.get_stage_options('trim_galore')
    	trim_args = 'trim_galore --paired --quality 20 --length 30 --stringency 5 --Phred33 ' \
	            '{input_files[0]} {input_files[1]} ' \
                    '--output_dir {output_paired_files[0]} {output_paired_files[1]} --clip_R1 25 --clip_R2 25 ' \
	            .format( read_in_R1 = input_files[0], read_in_R2=input_files[1],
		             out_R1=output_paired_files[0], 
                             out_R2= output_paired_files[1])
				
	
	run_stage(self.state, 'trim_galore', trim_args.format(**local()))

  	with logger_mutex:
       		logger.debug("Hooray trim_fastq worked")



    def align_bwa(self, inputs, bam_out, sample_id):
        '''Align the paired end fastq files to the reference genome using bwa'''
        fastq_read1_in, fastq_read2_in = inputs
        cores = self.get_stage_options('align_bwa', 'cores')
        read_group = '"@RG\tID:{sample}\tSM:{sample}\tPL:Illumina"'.format(sample=sample_id)
        command = 'bwa mem -t {cores} -R {read_group} {reference} {fastq_read1} {fastq_read2} ' \
                  '| samtools view -bS -h -o {bam} -' \
                  .format(cores=cores,
                      read_group=read_group,
                      fastq_read1=fastq_read1_in,
                      fastq_read2=fastq_read2_in,
                      reference=self.reference,
                      bam=bam_out)
        run_stage(self.state, 'align_bwa', command)


    def sort_bam_picard(self, bam_in, sorted_bam_out):
        '''Sort the BAM file using Picard'''
        picard_args = 'SortSam INPUT={bam_in} OUTPUT={sorted_bam_out} ' \
                      'VALIDATION_STRINGENCY=LENIENT SORT_ORDER=coordinate ' \
                      'CREATE_INDEX=True'.format(
                          bam_in=bam_in, sorted_bam_out=sorted_bam_out)
        self.run_picard('sort_bam_picard', picard_args)

