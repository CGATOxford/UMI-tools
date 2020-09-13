import numpy as np
cpdef int edit_distance(a, b):
    cdef char * aa = a
    cdef char * bb = b

    cdef int k, l, c

    c = 0

    la = len(a)
    # we only want to define hamming distances between barcodes with the same length 
    lb = len(b)
    if la != lb:
        return np.Inf
	    
    for k from 0 <= k < la:
        if aa[k] != bb[k]:
            c += 1
    return c
