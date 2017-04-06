'''
Build the pipeline workflow by plumbing the stages together.
'''

from ruffus import Pipeline, suffix, formatter, add_inputs, output_from
from stages import Stages
import ruffus
import sys
import functools
import decorator
import sys, os


#Using functools @wraps
def time_job(stream=sys.stdout):
    def actual_time_job(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            return time_func_call(func, stream, *args, **kwargs)
        return wrapper
    return actual_time_job

def time_job(stream=sys.stdout):
    def time_job(func, *args, **kwargs):
        return time_func_call(func, stream, *args, **kwargs)
    return decorator.decorator(time_job)

##By hand, using a callable object

class time_job(object):
    def __init__(self, stream=sys.stdout):
        self.stream = stream
    def __call__(self, func):
        def inner(*args, **kwargs):
            return time_func_call(func, self.stream, *args, **kwargs)
        # remember to forward __name__
        inner.__name__ = func.__name__
        inner.__module__ = func.__module__
        inner.__doc__ = func.__doc__
        if hasattr(func, "pipeline_task"):
            inner.pipeline_task = func.pipeline_task
        return inner

def time_func_call(func, stream, *args, **kwargs):
    """prints elapsed time to standard out, or any other file-like object with a .write() method. """
    start = time.time()
   # Run the decorated function.
    ret = func(*args, **kwargs)
   # Stop the timer.
    end = time.time()
    elapsed = end - start
    stream.write("{} took {} seconds\n".format(func.__name__, elapsed))
    return ret

def make_pipeline(state):
    '''Build the pipeline by constructing stages and connecting them together'''
    # Build an empty pipeline
    pipeline = Pipeline(name='complexo')
    # Get a list of paths to all the FASTQ files
    fastq_files = state.config.get_option('fastqs')
    # Find the path to the reference genome
    reference_file = state.config.get_option('ref_grch37')
    # Stages are dependent on the state
    stages = Stages(state)

    # The original reference file
    # This is a dummy stage. It is useful because it makes a node in the
    # pipeline graph, and gives the pipeline an obvious starting point.
    pipeline.originate(
         task_func=stages.original_reference,
         name='original_reference',
         output=reference_file)
 
    # The original FASTQ files
    # This is a dummy stage. It is useful because it makes a node in the
    # pipeline graph, and gives the pipeline an obvious starting point.
    
    pipeline.originate(
        task_func=stages.original_fastqs,
        name='original_fastqs',
        output=fastq_files)
    '''
    #Group together file pairs
    pipeline.collate(
	task_func=stages.original_fastqs,
        name='original_fastqs',
        # match file name up to the "R1.fastq.gz"
        filter=formatter(".+/(?P<sample>[a-zA-Z0-9]+)_R1.fastq.gz"),
        # Create output parameter supplied to next task
        output=["{path[0]}/{sample[0]}_R1.fastq.gz",  # paired file 1
        "{path[0]}/{sample[0]}_R2.fastq.gz"]) # paired file 2
    '''        
   
    # Index the reference using BWA 
    pipeline.transform(
        task_func=stages.index_reference_bwa,
        name='index_reference_bwa',
        input=output_from('original_reference'),
        filter=suffix('.fasta'),
        output=['.fasta.amb', '.fasta.ann', '.fasta.pac', '.fasta.sa', '.fasta.bwt'])
     
    # Index the reference using Novoindex 
    pipeline.transform(
        task_func=stages.index_reference_novoindex,
        name='index_reference_novoindex',
        input=output_from('original_reference'),
        filter=suffix('.fasta'),
        output='.nix')
    
    # Index the reference using samtools 
    pipeline.transform(
        task_func=stages.index_reference_samtools,
        name='index_reference_samtools',
        input=output_from('original_reference'),
        filter=suffix('.fasta'),
        output='.fasta.fai')
    
    # Trim paired end reads in FASTQ to the reference producing a BAM file
    pipeline.transform(
        task_func=stages.trim_galore,
        name='trim_galore',
        input=output_from('original_fastqs'),
	filter=formatter('.+/(?P<sample>[a-zA-Z0-9]+)_R1.fastq.gz'),
	add_inputs=add_inputs('{path[0]}/{sample[0]}_R2.fastq.gz'),
        #extras=['{sample[0]}'],
	#_R1_val_1.fq.gz
        output=r"'{path[0]}'/.fq.gz"
    )
             
  # Align paired end reads in FASTQ to the reference producing a BAM file
    
    (pipeline.transform(
        task_func=stages.novoalign,
        name='novoalign',
        input=output_from('trim_galore'),
        #filter=suffix('.fq.gz'),
        filter=formatter('.+/(?P<sample>[a-zA-Z0-9]+)_R1_val_1.fq.gz'),
        add_inputs=add_inputs(['{path[0]}/{sample[0]}_R2_val_2.fq.gz', reference_file]),
	   # sample specific configuration options
        extras=['{sample[0]}'],
        output="{path[0]}/{sample[0]}.bam")
	.follows('index_reference_novoindex'))
    ''' 
    # Sort the BAM file using Picard 
    pipeline.transform(
        task_func=stages.sort_bam_picard,
        name='sort_bam_picard',
        input=output_from('novoalign'),
        filter=suffix('.bam'),
        output='.sort.bam')


    # Mark duplicates in the BAM file using Picard 
    pipeline.transform(
        task_func=stages.mark_duplicates_picard,
        name='mark_duplicates_picard',
        input=output_from('sort_bam_picard'),
        filter=suffix('.sort.bam'),
        # xxx should make metricsup an extra output?
        output=['.sort.dedup.bam', '.metricsdup'])

    # Generate chromosome intervals using GATK 
    pipeline.transform(
        task_func=stages.chrom_intervals_gatk,
        name='chrom_intervals_gatk',
        input=output_from('mark_duplicates_picard'),
        filter=suffix('.sort.dedup.bam'),
        output='.chr.intervals')

    # Local realignment using GATK 
    (pipeline.transform(
        task_func=stages.local_realignment_gatk,
        name='local_realignment_gatk',
        input=output_from('chrom_intervals_gatk'),
        filter=formatter('.+/(?P<sample>[a-zA-Z0-9]+).chr.intervals'),
        add_inputs=add_inputs('{path[0]}/{sample[0]}.sort.dedup.bam'),
        output='{path[0]}/{sample[0]}.sort.dedup.realn.bam')
        .follows('mark_duplicates_picard'))

    # Base recalibration using GATK 
    pipeline.transform(
        task_func=stages.base_recalibration_gatk,
        name='base_recalibration_gatk',
        input=output_from('local_realignment_gatk'),
        filter=suffix('.sort.dedup.realn.bam'),
        output=['.recal_data.csv', '.count_cov.log'])

    # Print reads using GATK 
    (pipeline.transform(
        task_func=stages.print_reads_gatk,
        name='print_reads_gatk',
        input=output_from('base_recalibration_gatk'),
        filter=formatter('.+/(?P<sample>[a-zA-Z0-9]+).recal_data.csv'),
        add_inputs=add_inputs('{path[0]}/{sample[0]}.sort.dedup.realn.bam'),
        output='{path[0]}/{sample[0]}.sort.dedup.realn.recal.bam')
        .follows('local_realignment_gatk'))

    # Call variants using GATK 
    pipeline.transform(
        task_func=stages.call_variants_gatk,
        name='call_variants_gatk',
        input=output_from('print_reads_gatk'),
        filter=suffix('.sort.dedup.realn.recal.bam'),
        output='.raw.snps.indels.g.vcf')

    # Combine G.VCF files for all samples using GATK
    pipeline.merge(
        task_func=stages.combine_gvcf_gatk,
        name='combine_gvcf_gatk',
        input=output_from('call_variants_gatk'),
        output='PCExomes.mergegvcf.vcf')

    # SNP recalibration using GATK
    pipeline.transform(
        task_func=stages.snp_var_recalibrate_gatk,
        name='snp_var_recalibrate_gatk',
        input=output_from('call_variants_gatk'),
        filter=suffix('.raw.snps.indels.g.vcf'),
        output=['.snp_recal', '.snp_tranches', '.snp_plots.R'])

    # INDEL recalibration using GATK
    pipeline.transform(
        task_func=stages.indel_var_recalibrate_gatk,
        name='indel_var_recalibrate_gatk',
        input=output_from('call_variants_gatk'),
        filter=suffix('.raw.snps.indels.g.vcf'),
        output=['.indel_recal', '.indel_tranches', '.indel_plots.R'])

    # Apply SNP recalibration using GATK  
    (pipeline.transform(
        task_func=stages.apply_snp_var_recalibrate_gatk,
        name='apply_snp_var_recalibrate_gatk',
        input=output_from('call_variants_gatk'),
        filter=suffix('.raw.snps.indels.g.vcf'),
        add_inputs=add_inputs(['PCExomes.snp_recal', 'PCExomes.snp_tranches']),
        output='.recal_SNP.vcf')
        .follows('snp_var_recalibrate_gatk'))

    # Apply INDEL recalibration using GATK  
    (pipeline.transform(
        task_func=stages.apply_indel_var_recalibrate_gatk,
        name='apply_indel_var_recalibrate_gatk',
        input=output_from('call_variants_gatk'),
        filter=suffix('.raw.snps.indels.g.vcf'),
        add_inputs=add_inputs(['PCExomes.indel_recal', 'PCExomes.indel_tranches']),
        output='.recal_INDEL.vcf')
        .follows('indel_var_recalibrate_gatk'))

    # Combine variants using GATK  
    (pipeline.transform(
        task_func=stages.combine_variants_gatk,
        name='combine_variants_gatk',
        input=output_from('apply_snp_var_recalibrate_gatk'),
        filter=suffix('.recal_SNP.vcf'),
        add_inputs=add_inputs(['PCExomes.recal_INDEL.vcf']),
        output='.combined.vcf')
        .follows('apply_indel_var_recalibrate_gatk'))

    # Select variants using GATK 
    pipeline.transform(
        task_func=stages.select_variants_gatk,
        name='select_variants_gatk',
        input=output_from('combine_variants_gatk'),
        filter=suffix('.combined.vcf'),
        output='.selected.vcf')
    '''
    return pipeline
