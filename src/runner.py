'''
Run a pipeline stage either locally or on the cluster,
taking into account the configuration settings,
and command line options of the pipeline.
'''

from ruffus.drmaa_wrapper import run_job, error_drmaa_job

#Install drmaa  and set drmaa to virtual environment

#export SGE_ROOT="/var/lib/gridengine"

# sun grid engine memory is requested in MB, but the config file specifies in GB
MEGABYTES_IN_GIGABYTE = 1024

'''
SGE options:

-account=name          Charge job to specified accounts
-mem=MB                Minimum amount of real memory
-mem-per-cpu=MB        Maximum amount of real memory per allocated cpu
                       required by a job
-nodes=N               Number of nodes on which to run (N = min[-max])
-pe -serial=n            Number of tasks to invoke on each node
-M a	               Notify user by email when certain event types
                       occur. Valid type values are BEGIN, END, FAIL,
                       REQUEUE, and ALL (any state change)
'''

def run_stage(state, stage, command):
    '''Run a pipeline stage, either locally or on the cluster'''

    # Grab the configuration options for this stage
    config = state.config
    modules = config.get_stage_option(stage, 'modules')
    mem = config.get_stage_option(stage, 'mem') * MEGABYTES_IN_GIGABYTE 
    account = config.get_stage_option(stage, 'account')
    queue = config.get_stage_option(stage, 'queue')
    walltime = config.get_stage_option(stage, 'walltime')
    run_local = config.get_stage_option(stage, 'local')
    cores = config.get_stage_option(stage, 'cores')
    pipeline_id = config.get_option('pipeline_id')
    job_name = pipeline_id + '_' + stage

    # Generate a "module load" command for each required module
    
    module_loads = '. /etc/profile.d/module.sh\n' + '\n'.join(['module load ' + module for module in modules])
    cluster_command = '\n'.join([module_loads, command])

    #Specify job-specific options for SGE
    job_options = '-q {queue} -pe serial {cores} -l docker=1 -l h_vmem={mem}M  -M {account}'.format( cores=cores, queue=queue, mem=mem, account=account) 


    # Log a message about the job we are about to run
    log_messages = ['Running stage: {}'.format(stage),
                    'Command: {}'.format(command)]
    if not run_local:
        log_messages.append('Job options: {}'.format(job_options))
    state.logger.info('\n'.join(log_messages))

    # Run the job, capturing stdout and stderr
    stdout_res, stderr_res = None, None
    try:
        stdout_res, stderr_res = \
            run_job(cmd_str=cluster_command,
                job_name = job_name,
                logger = state.logger.proxy,
                drmaa_session = state.drmaa_session,
                # Determines whether to run the command on the local
                # machine or run it on the cluster
                run_locally = run_local,
                # Keep a copy of the job script for diagnostic purposes
                retain_job_scripts = True,
                #retain_stdout = True,
                #retain_stderr = True,
                job_script_directory = state.options.jobscripts, 
                job_other_options = job_options)
    except error_drmaa_job as err:
        raise Exception("\n".join(map(str, ["Failed to run:", command, err, stdout_res, stderr_res])))
