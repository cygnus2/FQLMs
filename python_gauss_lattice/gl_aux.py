""" ----------------------------------------------------------------------------

    gl_aux.py - LR, March 2020

    Holds some auxiliary stuff.

---------------------------------------------------------------------------- """
import time


def winding_tag(ws, labels=['x', 'y', 'z']):
    """ Returns the naming convention of the winding datasets.
    """
    wtag = ''
    for k in range(len(ws)):
        wtag += 'w{:s}_{:d}-'.format(labels[k], ws[k])
    return wtag[:-1]


def file_tag(L, filetype='hdf5'):
    """ Returns the naming convention of the winding datasets.
    """
    ftag = ''
    for k in range(len(L)):
        ftag += '{:d}x'.format(L[k])
    return 'winding_states_' + ftag[:-1] + '.' + filetype



def timeit(method):
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        if 'log_time' in kw:
            name = kw.get('log_name', method.__name__.upper())
            kw['log_time'][name] = int((te - ts) * 1000)
        else:
            print ('%r  %2.2f ms' % (method.__name__, (te - ts) * 1000))
        return result
    return timed


def print_2D_state(state, L):
    """ Dumps a 2D lattice.
    """
    sx, sy = 2, 1
    blx = 2*(sx+1)

    bstr = bin_str(state, L=L**2*2)[::-1]
    bstr = bstr.replace('0', '○')
    bstr = bstr.replace('1', '●')

    # Loop through rows.
    lines = []
    for r in range(L):
        lx = ''
        for c in range(L):
            i = 2*r*L+c*2
            lx += "-"+"-"*sx + bstr[i] + "-"*sx
        lines.append(lx)

        # ---

        for k in range(sy):
            lines.append(("|"+(" "*blx)[:-1])*L)

        ly = ''
        for c in range(L):
            j = 2*r*L+c*2 + 1
            ly += bstr[j]+(" "*blx)[:-1]
        lines.append(ly)

        for k in range(sy):
            lines.append(("|"+(" "*blx)[:-1])*L)

    for line in lines[::-1]:
        print(line)
