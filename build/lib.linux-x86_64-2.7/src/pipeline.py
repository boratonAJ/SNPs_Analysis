'''
Build the pipeline workflow by plumbing the stages together.
'''

from ruffus import Pipeline, suffix, formatter, add_inputs, output_from
from stages import Stages
import ruffus
import sys

def make_pipeline(state):
    '''Build the pipeline by constructing stages and connecting them together'''
    # Build an empty pipeline
    pipeline = Pipeline(name='complexo')
    # Get a list of paths to all the FASTQ files
    fastq_files = state.config.get_option('fastqs')
    # Find the path to the reference genome
    reference_file = state.config.get_option('ref_grch37')
    index_file = state.config.get_option('index_file')
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
    
    #Create dictionary for the reference using picard_tool 
    pipeline.transform(
        task_func=stages.createSequenceDictionary,
        name='createSequenceDictionary',
        input=output_from('original_reference'),
        filter=suffix('.fasta'),
        output='.dict')

    #Trim for data quality
    pipeline.transform(
	task_func= stages.trim_fastq,
	name='trim_fastq',
	# match file name up to the "R1.fastq.gz"
	input=output_from('original_fastqs'),
        filter=formatter('.+/(?P<Test>[a-zA-Z0-9]+)_R1.fastq.gz'), #'{path[0]}/{sample[0]}_R2.fastq.gz'),
	add_inputs=add_inputs('{path[0]}/{Test[0]}_R2.fastq.gz', '{path[0]}/'),
	# Extra parameters for our own convenience and use
	extras=['{Test[0]}'],
        #extra=['"{path[0]}/{1[0]}unpaired.R1.fastq.gz","{path[0]}/{1[0]}unpaired.R2.fastq.gz"'], # unpaired file 2
        # Create output parameter supplied to next task
        output=['{path[0]}/{1[0]}_1P_R1.fq.gz', '{path[0]}/{1[0]}_2P_R2.fq.gz',# forward & reverse reads 
	'{path[0]}/{1[0]}_1UP_R1.fq.gz', '{path[0]}/{1[0]}_2UP_R2.fq.gz'] # unpaired file 1&2
        )
      
     # Novo_Align paired end reads in FASTQ to the reference producing a BAM file
    (pipeline.transform(
        task_func=stages.novoalign,
        name='novoalign',
        input=output_from('trim_fastq'),
        #filter=suffix(".fq.gz"),
	filter = formatter('.+/(?P<Test>[a-zA-Z0-9]+)_[0-9]P_R[0-9].fq.gz'),
        #filter=formatter(['.+/(?P<sample>\w+)\_R2_val_2.fq.gz', reference_file]),
        add_inputs=add_inputs([index_file]),
	   # sample specific configuration options
        extras=['{Test[0]}'],
        output="{path[0]}/{Test[0]}.bam")
	.follows('index_reference_novoindex')
	.follows('index_reference_samtools'))
    
     
    # Sort the BAM file using novosort 
    pipeline.transform(
        task_func=stages.novosort_bam,
        name='novosort_bam',
        input=output_from('novoalign'),
        filter=suffix('.bam'),
        output='.sort.bam')
    
    # Align paired end reads in FASTQ to the reference producing a BAM file
    (pipeline.transform(
        task_func=stages.align_bwa,
        name='align_bwa',
        input=output_from('original_fastqs'),
        filter=formatter('.+/(?P<sample>[a-zA-Z0-9]+)_R1.fastq.gz'),
        add_inputs=add_inputs('{path[0]}/{sample[0]}_R2.fastq.gz'),
        extras=['{sample[0]}'],
        output='{path[0]}/{sample[0]}.bwa.bam')
        .follows('index_reference_bwa')
        .follows('index_reference_samtools'))

    # Sort the BAM file using Picard 
    pipeline.transform(
        task_func=stages.sort_bam_picard,
        name='sort_bam_picard',
        input=output_from('align_bwa'),
        filter=suffix('.bam'),
        output='.sort.bwa.bam')    

    # Mark duplicates in the BAM file using Picard 
    pipeline.transform(
        task_func=stages.mark_duplicates_picard,
        name='mark_duplicates_picard',
        input=output_from('novosort_bam'),
        filter=suffix('.sort.bam'),
        # xxx should make metricsup an extra output?
        output=['.sort.dedup.bam', '.metricsdup'])
    
 # Generate an alignment statistics report for the sorted bam file using picard tools
    pipeline.transform(
        task_func=stages.collectalignmentmetric,
        name='collectalignmentmetric',
        input=output_from('mark_duplicates_picard'),
        filter=suffix('.sort.dedup.bam'),
        #add_inputs=add_inputs([reference_file]),
        output='.stat.txt')

    pipeline.transform(
        task_func=stages.buildBamIndex,
        name='buildBamIndex',
        input=output_from('mark_duplicates_picard'),
        filter=suffix('.sort.dedup.bam'),
        # xxx should make metricsup an extra output?
        output='.bai')
    
    # local base realignertargetCreator using GATK 
    pipeline.transform(
        task_func=stages.RealignerTarget_gatk,
        name='RealignerTarget_gatk',
        input=output_from('mark_duplicates_picard'),
        filter=suffix('.sort.dedup.bam'),
        output='.chr.intervals')
     
    # Local base realignment using GATK 
    (pipeline.transform(
        task_func=stages.IndelRealigner_gatk,
        name='IndelRealigner_gatk',
        input=output_from('RealignerTarget_gatk'),
        filter=formatter('.+/(?P<Test>[a-zA-Z0-9]+).chr.intervals'),
        add_inputs=add_inputs('{path[0]}/{Test[0]}.sort.dedup.bam'),
        output='{path[0]}/{Test[0]}.sort.dedup.realn.bam')
        .follows('mark_duplicates_picard'))
    '''
    # local base recalibration using GATK 
    pipeline.transform(
        task_func=stages.basedRecalibrate_gatk,
        name='basedRecalibrate_gatk',
        input=output_from('IndelRealigner_gatk'),
        filter=suffix('.sort.dedup.realn.bam'),
        output=['.recal_data.csv', '.count_cov.log'])

    # Print reads using GATK 
    (pipeline.transform(
        task_func=stages.viewReads_gatk,
        name='viewReads_gatk',
        input=output_from('basedRecalibrate_gatk'),
        filter=formatter('.+/(?P<Test>[a-zA-Z0-9]+).recal_data.csv'),
        add_inputs=add_inputs('{path[0]}/{Test[0]}.sort.dedup.realn.bam'),
        output='{path[0]}/{Test[0]}.sort.dedup.realn.recal.bam')
        .follows('IndelRealigner_gatk'))
    '''
    # Call variants using GATK 
    pipeline.transform(
        task_func=stages.call_variants_gatk,
        name='call_variants_gatk',
        input=output_from('IndelRealigner_gatk'),
        filter=suffix('.sort.dedup.realn.bam'),
        output='.raw.snps.indels.g.vcf')
    
    pipeline.transform(
        task_func=stages.call_variants_samtools,
        name='call_variants_samtools',
        input=output_from('IndelRealigner_gatk'),
        filter=suffix('.sort.dedup.realn.bam'),
        output='.raw.sam_bwa.snps.indels.g.vcf')

    pipeline.transform(
        task_func=stages.variant_samtools_filter,
        name='variant_samtools_filter',
        input=output_from('call_variants_samtools'),
        filter=suffix('.raw.sam_bwa.snps.indels.g.vcf'),
        output= '_filtered.vcf')
   
    pipeline.transform(
        task_func=stages.freebayes_call,
        name='freebayes_call',
        input=output_from('IndelRealigner_gatk'),
        filter=suffix('.sort.dedup.realn.bam'),
        output='.raw.freebayes.snps.indels.g.vcf')

    # Call snpeff using SNPEFF Variants 
    pipeline.transform(
        task_func=stages.annotating_vcf_SNPEFF,
        name='annotating_vcf_SNPEFF',
        input=output_from('call_variants_gatk'),
        filter=suffix('.raw.snps.indels.g.vcf'),
        output='.vcf')


    return pipeline
