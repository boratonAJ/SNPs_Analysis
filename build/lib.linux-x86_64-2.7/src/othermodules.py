#This function use a function from another module as a decorator, 
#but I need it to manipulate the current module's global namespace.

import sys
def decorator(cls):
    mod = sys.modules[cls.__module__]
    mod.root = cls
