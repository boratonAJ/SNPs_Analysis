#!usr/bin/en python

#What would @time_job look like?
#Using functools @wraps


from ruffus import *
import sys
import time
import decorator

def time_job(stream=sys.stdout):
    def time_job(func, *args, **kwargs):
        return time_func_call(func, stream, *args, **kwargs)
    return decorator.decorator(time_job)

'''Using Michele Simionatos decorator module'''

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
    """prints elapsed time to standard out, or any other file-like object with a .write() method.
    """
    start = time.time()
   # Run the decorated function.
    ret = func(*args, **kwargs)
   # Stop the timer.
    end = time.time()
    elapsed = end - start
    stream.write("{} took {} seconds\n".format(func.__name__, elapsed))
    return ret



@time_job(sys.stderr)
def first_task():
    print "First task"


@follows(first_task)
@time_job(sys.stderr)
def second_task():
    print "Second task"


@follows(second_task)
@time_job(sys.stderr)
def final_task():
    print "Final task"

pipeline_run()



