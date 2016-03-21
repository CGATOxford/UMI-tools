cpdef int edit_distance(a, b):
    cdef char * aa = a
    cdef char * bb = b

    cdef int k, l, c

    c = 0

    l = len(a)
    for k from 0 <= k < l:
        if aa[k] != bb[k]:
            c += 1
    return c
