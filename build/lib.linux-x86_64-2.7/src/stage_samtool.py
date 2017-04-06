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
import sys
import glob


PICARD_JAR='$PICARD_HOME/picard.jar'
GATK_JAR='$GATK_HOME/GenomeAnalysisTK.jar'
TRIM_JAR='/cip0/software/java/Trimmomatic-0.33/trimmomatic-0.33.jar'
SNPEFF_JAR='/usr/people/ajayi/test/complexo_pipeline/snpEff/snpEff.jar'

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

    def run_picard(self, stage, args):
        mem = int(self.state.config.get_stage_options(stage, 'mem'))
        return run_java(self.state, stage, PICARD_JAR, mem, args)

    def run_gatk(self, stage, args):
        mem = int(self.state.config.get_stage_options(stage, 'mem'))
        return run_java(self.state, stage, GATK_JAR, mem, args)
    
    def run_trimomatic(self, stage, args):
        mem = int(self.state.config.get_stage_options(stage, 'mem'))
	return run_java(self.state, stage, TRIM_JAR, mem, args)

    def run_snpeff(self, stage, args):
        mem = int(self.state.config.get_stage_options(stage, 'mem'))
        return run_java(self.state, stage, SNPEFF_JAR, mem, args)

    
    def get_stage_options(self, stage, *options):
        return self.state.config.get_stage_options(stage, *options)

    def get_options(self, *options):
        return self.state.config.get_options(*options)

    def original_reference(self, output):
        '''Original reference file'''
        pass

    def original_fastqs(self, output):
        '''Original fastq files'''
        pass
    
    def index_reference_bwa(self, reference_in, index_file_out):
        '''Index the reference genome using BWA'''
        command = "bwa index -a bwtsw {ref}".format(ref=reference_in)
        run_stage(self.state, 'index_reference_bwa', command)
    
    def index_reference_samtools(self, reference_in, index_file_out):
        '''Index the reference genome using samtools'''
        command = 'samtools faidx {ref}'.format(ref=reference_in)
        run_stage(self.state, 'index_reference_samtools', command)
    
    def createSequenceDictionary(self, reference_in, outputs):
	'''Creating reference dictionary using Picard'''
        dict_out = outputs
        picard_args = 'CreateSequenceDictionary REFERENCE={ref} OUTPUT={dict_out} ' \
                      'VALIDATION_STRINGENCY=LENIENT ' \
                      'CREATE_INDEX=True'.format(ref=reference_in, dict_out=dict_out,
                          )
        self.run_picard('createSequenceDictionary', picard_args)

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

    def mark_duplicates_picard(self, bam_in, outputs):
        '''Mark duplicate reads using Picard'''
        dedup_bam_out, metrics_out = outputs
        picard_args = 'MarkDuplicates INPUT={bam_in} OUTPUT={dedup_bam_out} ' \
                      'METRICS_FILE={metrics_out} VALIDATION_STRINGENCY=LENIENT ' \
                      'MAX_RECORDS_IN_RAM=5000000 ASSUME_SORTED=True ' \
                      'COMPRESSION_LEVEL=5 CREATE_INDEX=True'.format(bam_in=bam_in, dedup_bam_out=dedup_bam_out,
                          metrics_out=metrics_out)
        self.run_picard('mark_duplicates_picard', picard_args)

    def collectalignmentmetric(self, inputs, output):
        '''Collect stats from bam file using Picard'''
	#print ">>>>>>>>CollectMetrisInput<<<<<<<<<<", inputs[1][0]
	bam_in, _metrics_dup = inputs
        picard_args = 'CollectAlignmentSummaryMetrics REFERENCE_SEQUENCE={reference} INPUT={bam_input} OUTPUT={bam_output}' \
                      .format(bam_input= bam_in, reference=self.reference, bam_output = output)
        self.run_picard('collectalignmentmetric', picard_args)
    
    #Index the Bam files
    def buildBamIndex(self, inputs, outputs):
        '''Mark duplicate reads using Picard'''
        dedup_bam_in, metrics_in = inputs
        picard_args = 'BuildBamIndex INPUT={bam_in} OUTPUT={bai_out} ' \
                      .format(bam_in=dedup_bam_in, bai_out=outputs)
        self.run_picard('buildBamIndex', picard_args)
    
    # RealignerTargetCreator-Based Realignment:
    def RealignerTarget_gatk(self, inputs, intervals_out):
        '''Generate chromosome intervals using GATK'''
        bam_in, _metrics_dup = inputs
        cores = self.get_stage_options('RealignerTarget_gatk', 'cores')
        gatk_args = '-T RealignerTargetCreator -R {reference} -I {bam} ' \
                    '--num_threads {threads} ' \
                    '-o {out}'.format(reference=self.reference, bam=bam_in,
                            threads=cores, out=intervals_out)
        self.run_gatk('RealignerTarget_gatk', gatk_args)

    #Indel-Based Realigner
    def IndelRealigner_gatk(self, inputs, realigned_bam_out):
        '''Local realign reads using GATK''REFERENCE'''
        target_intervals_in, bam_in = inputs
	cores = self.get_stage_options('IndelRealigner_gatk', 'cores')
        gatk_args = "-T IndelRealigner -R {reference} -I {bam} " \
                    "-targetIntervals {target_intervals} " \
                    "-o {out}".format(reference=self.reference, bam=bam_in,
                            target_intervals=target_intervals_in,
                            out=realigned_bam_out)
        self.run_gatk('IndelRealigner_gatk', gatk_args)
   
    # XXX I'm not sure that --num_cpu_threads_per_data_thread has any benefit here
    def basedRecalibrate_gatk(self, bam_in, outputs):
        '''Base recalibration using GATK'''
        csv_out, log_out = outputs
        gatk_args = "-T BaseRecalibrator -R {reference} -I {bam} " \
                    "--num_cpu_threads_per_data_thread 4 " \
                    "-log {log} -o {out}".format(reference=self.reference, bam=bam_in,
                            log=log_out, out=csv_out)
        self.run_gatk('basedRecalibrate_gatk', gatk_args)


    def viewReads_gatk(self, inputs, bam_out):
        '''Print reads using GATK'''
        [csv_in, _log], bam_in = inputs
        gatk_args = "-T PrintReads -R {reference} -I {bam} --BQSR {recal_csv} " \
                    "-o {out} --num_cpu_threads_per_data_thread 4".format(reference=self.reference,
                            bam=bam_in, recal_csv=csv_in, out=bam_out)
        self.run_gatk('viewReads_gatk', gatk_args)


    def call_variants_gatk(self, bam_in, vcf_out):
        '''Call variants using GATK'''
        gatk_args = "-T HaplotypeCaller -R {reference} --min_base_quality_score 20 " \
                    "--variant_index_parameter 128000 --emitRefConfidence GVCF " \
                    "--standard_min_confidence_threshold_for_calling 30.0 " \
                    "--num_cpu_threads_per_data_thread 8 " \
                    "--variant_index_type LINEAR " \
                    "--standard_min_confidence_threshold_for_emitting 30.0 " \
                    "-I {bam} -o {out}".format(reference=self.reference,
                            bam=bam_in, out=vcf_out)
        self.run_gatk('call_variants_gatk', gatk_args)



