'''
Build the pipeline workflow by plumbing the stages together.
'''

from ruffus import Pipeline, suffix, formatter, add_inputs, output_from
from stages import Stages


def make_pipeline(state):
    '''Build the pipeline by constructing stages and connecting them together'''
    # Build an empty pipeline
    pipeline = Pipeline(name='complexo')
    # Get a list of paths to all the FASTQ files
    fastq_files = state.config.get_option('fastqs')
    # Stages are dependent on the state
    stages = Stages(state)
	
    # The original FASTQ files
    # This is a dummy stage. It is useful because it makes a node in the
    # pipeline graph, and gives the pipeline an obvious starting point.
    pipeline.originate(
        task_func=stages.original_fastqs,
        name='original_fastqs',
        output= fastq_files)


    # Trimming paired end reads in FASTQ to get quality file
    pipeline.transform(
        task_func=stages.trim_galore,
        name='trim_galore',
        input=output_from('original_fastqs'),
        # characters.sample1_R2.fastq.gz
        filter=formatter('.+/(?P<sample>[a-zA-Z0-9]+)_R1.fastq.gz'),
        # Add one more inputs to the stage:
        #    1. The corresponding R2 FASTQ file
        #add_inputs=add_inputs('{path[0]}/{sample[0]}_R2.fastq.gz'),
	#extras=['{sample[0]}'],
       #output= "{path[0]}")
	output=["{path[0]}/{1[0]}_R1_val_1.fq.gz",
		"{path[1]}/{1[0]}_R2_val_2.fq.gz"])


    return pipeline
