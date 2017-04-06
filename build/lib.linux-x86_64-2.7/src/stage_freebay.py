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
    
    def index_reference_novoindex(self, reference_in, index_file_out):
        '''Index the reference genome using novoindex'''
        command = "novoindex {out} {ref}".format(ref=reference_in, out=index_file_out)
        run_stage(self.state, 'index_reference_novoindex', command)
    
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

    def trim_fastq(self, input_files, output_paired_files, discarded_unpaired_files):
	#start = time.time()
	if len(input_files) != 3:
        	raise Exception("inputs:", input_files)
	input_files = sorted((glob.glob(input_files[2] + '*.fastq.gz')))
	fq_read1_in = input_files[0] 
	fq_read2_in= input_files[1]
	out_fasta_PE_1, out_fasta_PE_2 = output_paired_files[:2]
	out_fasta_U_1, out_fasta_U_2 = output_paired_files[2:]
	cores = self.state.config.get_stage_option('trim_fastq', 'cores') 
        trim_arg ="PE -threads {cores} -phred33 {input_files1} {input_files2} " \
           	  " {output_paired_files1} {discarded_unpaired_files1} " \
           	  " {output_paired_files2} {discarded_unpaired_files2} " \
           	  " LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 " \
	   	  .format (cores=cores,
		      input_files1= fq_read1_in,
		      input_files2= fq_read2_in,
		      output_paired_files1 = out_fasta_PE_1,
		      output_paired_files2 = out_fasta_PE_2,
		      discarded_unpaired_files1 = out_fasta_U_1,
		      discarded_unpaired_files2 = out_fasta_U_2	
			)
	self.run_trimomatic('trim_fastq', trim_arg)

    #Align Illumina paired-end reads to a reference genome. Specify the expected size distriibution 
    def novoalign(self, inputs, bam_out, sample):
        '''Align the paired end fastq files to the reference genome using novoalign '''
	fastq_read1_in, fastq_read2_in = inputs[0][:2]
	reference_in = inputs[1][0]
        # Get the number of cores to request for the job, this translates into the
        # number of threads to give to bwa's -t option
        cores = self.state.config.get_stage_option('novoalign', 'cores')
        read_group = '"@RG\tID:{sample}\tSM:{sample}\tPL:Illumina"'.format(sample=sample)
        # Run novoalign and pipe the output through samtools view to generate a BAM file
        command = 'novoalign -c {cores} -k -d {reference} -f {fastq_read1} {fastq_read2} ' \
		  '-i 250,100 -o SAM {read_group} ' \
                  '| samtools view -bS -o {bam} -' \
                  .format(cores=cores,
                      read_group=read_group,
                      fastq_read1=fastq_read1_in,
                      fastq_read2=fastq_read2_in,
                      reference=reference_in,
                      bam=bam_out)
        run_stage(self.state, 'novoalign', command)
	
    #novosort
    def novosort_bam(self, bam_in, sorted_bam_out):
        '''Sort the BAM file using novosort'''
	cores = self.state.config.get_stage_option('novoalign', 'cores')
	#mem = int(self.state.config.get_stage_options(stage, 'mem'))
        command = 'novosort -c {cores} -m 10G -s -f {bam_in} -o {sorted_bam_out} ' \
     		  .format( cores = cores, bam_in=bam_in, sorted_bam_out=sorted_bam_out)
        run_stage(self.state, 'novosort_bam', command)

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

    # Annotate .vcf files using SnpEff
    def annotating_vcf_SNPEFF(self, inputs, outputs):
        '''Annotating vcf files using output of variant calling'''
        vcf_in = inputs
	#cores = self.get_stage_options('annotating_vcf_SNPEFF', 'cores')
        snpeff_args = '{reference} {vcf_in} > {anno_out} ' \
                      .format(reference=self.reference, vcf_in=vcf_in, anno_out=outputs)
        self.run_snpeff('annotating_vcf_SNPEFF', snpeff_args)
        #run_stage(self.state, 'annotating_vcf_SNPEFF', snpeff_args)

    def combine_gvcf_gatk(self, vcf_files_in, vcf_out):
        '''Combine G.VCF files for all samples using GATK'''
        g_vcf_files = ' '.join(['--variant ' + vcf for vcf in vcf_files_in])
        gatk_args = "-T CombineGVCFs -R {reference} " \
                    "--disable_auto_index_creation_and_locking_when_reading_rods " \
                    "{g_vcf_files} -o {vcf_out}".format(reference=self.reference,
                        g_vcf_files=g_vcf_files, vcf_out=vcf_out)
        self.run_gatk('combine_gvcf_gatk', gatk_args)

    def snp_var_recalibrate_gatk(self, raw_snps_indels_g_vcf_in, outputs):
        '''SNP recalibration using GATK'''
        recal_snp_out, tranches_snp_out, snp_plots_r_out = outputs
        cores = self.get_stage_options('snp_var_recalibrate_gatk', 'cores')
        gatk_args = "-T VariantRecalibrator --disable_auto_index_creation_and_locking_when_reading_rods " \
                    "-R {reference} --minNumBadVariants 5000 --num_threads {cores} " \
                    "-resource:hapmap,known=false,training=true,truth=true,prior=15.0 {hapmap} " \
                    "-resource:omni,known=false,training=true,truth=false,prior=12.0 {one_k_g_snps} " \
                    "-resource:1000G,known=false,training=true,truth=false,prior=10.0 {one_k_g_highconf_snps} " \
		    "-resource:dbsnp,known=true,training=false,truth=false,prior=8.0 {dbsnp_grch37} " \
                    "-an ReadPosRankSum -an FS -an DP -an InbreedingCoeff " \
                    "-input {raw_snps_indels_g_vcf} --recal_file {recal_snp} --tranches_file {tranches_snp} " \
                    "-rscriptFile {snp_plots} -mode SNP".format(reference=self.reference, cores=cores, hapmap=self.hapmap, one_k_g_snps=self.one_k_g_snps, one_k_g_highconf_snps=self.one_k_g_highconf_snps, dbsnp_grch37=self.dbsnp_grch37, raw_snps_indels_g_vcf=raw_snps_indels_g_vcf_in, recal_snp=recal_snp_out, tranches_snp=tranches_snp_out, snp_plots=snp_plots_r_out)
        self.run_gatk('snp_var_recalibrate_gatk', gatk_args)


    def indel_var_recalibrate_gatk(self, raw_snps_indels_g_vcf_in, outputs):
        '''INDEL recalibration using GATK'''
        recal_indel_out, tranches_indel_out, indel_plots_r_out = outputs
        cores = self.get_stage_options('indel_var_recalibrate_gatk', 'cores')
        gatk_args = "-T VariantRecalibrator --disable_auto_index_creation_and_locking_when_reading_rods " \
                    "-R {reference} --minNumBadVariants 5000 --num_threads {cores} " \
                    "-resource:mills,known=false,training=true,truth=true,prior=12.0 {mills_grch37} " \
                    "-resource:1000G,known=false,training=true,truth=true,prior=10.0 {one_k_g_indels} " \
                    "-an ReadPosRankSum -an FS -an DP -an InbreedingCoeff " \
                    "-input {raw_snps_indels_g_vcf} -recalFile {recal_indel} " \
                    "-tranchesFile {tranches_indel} -rscriptFile {indel_plots} " \
                    " -mode INDEL".format(reference=self.reference,
                            cores=cores, mills_grch37=self.mills_grch37, one_k_g_indels=self.one_k_g_indels,
                            raw_snps_indels_g_vcf=raw_snps_indels_g_vcf_in, recal_indel=recal_indel_out,
                            tranches_indel=tranches_indel_out, indel_plots=indel_plots_r_out)
        self.run_gatk('indel_var_recalibrate_gatk', gatk_args)


    def apply_snp_var_recalibrate_gatk(self, inputs, vcf_out):
        '''Apply SNP recalibration using GATK'''
        genotype_vcf_in, [recal_snp, tranches_snp] = inputs 
        cores = self.get_stage_options('apply_snp_var_recalibrate_gatk', 'cores')
        gatk_args = "-T ApplyRecalibration --disable_auto_index_creation_and_locking_when_reading_rods " \
                    "-R {reference} --ts_filter_level 99.5 --excludeFiltered --num_threads {cores} " \
                    "-input {raw.snps.indels.g.vcf} -recalFile {recal_snp} -tranchesFile {tranches_snp} " \
                    "-mode SNP -o {vcf_out}".format(reference=self.reference,
                            cores=cores, raw_snps_indels_g_vcf=raw_snps_indels_g_vcf_in, recal_snp=recal_snp,
                            tranches_snp=tranches_snp, vcf_out=vcf_out)
        self.run_gatk('apply_snp_var_recalibrate_gatk', gatk_args)


    def apply_indel_var_recalibrate_gatk(self, inputs, vcf_out):
        '''Apply INDEL recalibration using GATK'''
        genotype_vcf_in, [recal_indel, tranches_indel] = inputs 
        cores = self.get_stage_options('apply_indel_var_recalibrate_gatk', 'cores')
        gatk_args = "-T ApplyRecalibration --disable_auto_index_creation_and_locking_when_reading_rods " \
                    "-R {reference} --ts_filter_level 99.0 --excludeFiltered --num_threads {cores} " \
                    "-input {genotype_vcf} -recalFile {recal_indel} -tranchesFile {tranches_indel} " \
                    "-mode INDEL -o {vcf_out}".format(reference=self.reference,
                            cores=cores, raw_snps_indels_g_vcf=raw_snps_indels_g_vcf_in, recal_indel=recal_indel,
                            tranches_indel=tranches_indel, vcf_out=vcf_out)
        self.run_gatk('apply_indel_var_recalibrate_gatk', gatk_args)


    def combine_variants_gatk(self, inputs, vcf_out):
        '''Combine variants using GATK'''
        recal_snp, [recal_indel] = inputs 
        cores = self.get_stage_options('combine_variants_gatk', 'cores')
        gatk_args = "-T CombineVariants -R {reference} --disable_auto_index_creation_and_locking_when_reading_rods " \
                    "--num_threads {cores} --genotypemergeoption UNSORTED --variant {recal_snp} " \
                    "--variant {recal_indel} -o {vcf_out}".format(reference=self.reference,
                            cores=cores, recal_snp=recal_snp, recal_indel=recal_indel,
                            vcf_out=vcf_out)
        self.run_gatk('combine_variants_gatk', gatk_args)


    def select_variants_gatk(self, combined_vcf, vcf_out):
        '''Select variants using GATK'''
        gatk_args = "-T SelectVariants -R {reference} --disable_auto_index_creation_and_locking_when_reading_rods " \
                    "--variant {combined_vcf} -select 'DP > 100' -o {vcf_out}".format(reference=self.reference,
                            combined_vcf=combined_vcf, vcf_out=vcf_out)
        self.run_gatk('select_variants_gatk', gatk_args)
