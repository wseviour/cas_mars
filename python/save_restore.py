"""
This module provides similar functionality to
IDL save/restore.

Adapted from:
http://idl2python.blogspot.co.uk/2010/10/save-and-restore-2.html

Author: William Seviour
Last update: 21/9/12
"""
def save(file,**kwargs):
    """
    Save the value of some data in a file.
    Usage: save('misdatos.pypic',a=a,b=b,test=test)
    """
    import pickle
    f=open(file,"wb")
    pickle.dump(kwargs,f,protocol=2)
    f.close

def restore(file):
    """
    Read data saved with save function.
    Usage: datos = restore('misdatos.pypic')
    """
    import pickle
    f=open(file,"rb")
    result = pickle.load(f)
    f.close
    return result
