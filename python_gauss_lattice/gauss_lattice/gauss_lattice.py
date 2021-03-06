""" ----------------------------------------------------------------------------

    gauss_lattice.py - LR, March 2020

    Implements a recursive algorithm to find states on a d-dimensional lattice
    that obey Gauss' Law.

---------------------------------------------------------------------------- """
import numpy as np
import math, os, datetime
import h5py as hdf
from itertools import product
from queue import Empty, Queue
from multiprocessing import Pool
from .aux_stuff import winding_tag, charge_tag
import sys

class GaussLattice(object):
    """ Represents a Gauss lattice in arbitrary dimension.

        To efficiently loop through the lattice states, it is convenient to pre-
        compute various things. The rough strategy is as follows:

            [Essential]
            1)  Compute the basis vertices (that obey GL) on a d-dimensional
                square lattice
            2)  Compute the lattice version of the basis to save effort later.
            3)  Pre-compute indicies/masks of states to check
            4)  Get the list of states that are able to rule out a configuration
                with a certain amount of information (checkable states with partial
                basis strings)

        Notes:
            -   tried a string representation, turns out this is ~10% slower than
                the index-based representation of the basis states.
    """

    def __init__(self, L, state_file=None, **kwargs):
        self.L = np.array(L)
        self.d = len(L)

        # Compute partial volume factor & find number of sites in the sublattice.
        self.S = [1]
        for l in L:
            self.S.append(l*self.S[-1])
        self.N_sublattice = self.S[-1] // 2
        self.bitlen = self.S[-1]*self.d

        # Binary string format for debugging purposes.
        self.binformat = '{:0'+str(self.N_sublattice*4)+'b}'

        # Make the charge background (for GL checks).
        self.static_charge_indicies = [
            [c-1 for c in clist] for clist in kwargs.get('static_charges', [[],[]])
        ] # conversion from 1-based notation.
        self.static_charges = self.make_charge_background(self.static_charge_indicies)
        self.has_charges = any(self.static_charges)
        print("static charge background: ", self.static_charges)

        # These are the positions of the vertices that belong to the sublattices.
        ind_bs, ind_gls = self.checkerboard_indices()
        assert len(ind_bs) == self.N_sublattice
        assert len(ind_gls) == self.N_sublattice

        # Construct the basis vertex for the lattice / charge configuration.
        vertex_base = self.find_vertex_base()
        self.spatial_base_elements = self.create_spatial_base_elements(ind_bs, self.static_charges, vertex_base)
        for i, bs in enumerate(ind_bs):
            if self.static_charges[bs] == 0:
                assert len(self.spatial_base_elements[i]) == math.factorial(2*self.d)/math.factorial(self.d)**2

        # Construct the basis vertices on the sublattice A. This reduces the cost
        # for later operation since no arithmetic must be done except for a single
        # lookup.
        self.lattice_base = self.construct_lattice_basis(ind_bs, self.static_charges, vertex_base)
        assert len(self.lattice_base) == len(ind_bs)

        # Pre-compute stuff to check the validity of Gauss Law.
        self.mpos, self.mneg = self.create_gls_masks(ind_gls)
        self.checkable_gls = self.find_checkable_gls(ind_bs, ind_gls)
        assert len(self.checkable_gls[-1]) == self.N_sublattice


        # Initialize winding number bins.
        self.use_winding = not self.has_charges
        if self.use_winding:
            self.winding_bins, self.winding_masks = self.prepare_winding_numbers()
        else:
            self.winding_bins = np.zeros(shape=(1,),dtype=np.int)
            if self.has_charges:
                self.ds_label = charge_tag(self.static_charge_indicies)
            else:
                self.ds_label = "all-ws"


        # ----------------------------------------------------------------------
        # Some settings for storage and I/O.

        # For winding number counts.
        self.filename ='winding_sectors{:d}.dat'
        self.out_dir = kwargs.get('basedir')
        if self.out_dir is None:
            self.out_dir = 'output/'
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)


        # "state_file" dictates the output behavior:
        #  - if set, states and winding numbers will be written to file
        #  - the file extension dictates the type:
        #       - .hdf5 for HDF5
        #       - .* for plain text file.
        self.write_states = bool(state_file)
        if self.write_states:
            # Writing to file should happen through a buffer, since otherwises
            # we'd need to do a costly file I/O operation after every state,
            # which greatly slows down the calculation. In order to circumvent
            # race conditions in the recursion, the buffer is realized as a queue,
            # which handles multi-producer data quite well (the order does not
            # really matter for us anyways). There's some price to pay in terms
            # of runtime, but it's generally not the bottleneck.
            buf_len = kwargs.get('buffer_length')
            if buf_len is None:
                buf_len = 1e6
            self.buffer = Queue(maxsize=buf_len)

            # Initialize the file.
            append = kwargs.get('append_states', False)
            print(append)
            self._init_file(state_file, append=append)


    @staticmethod
    def _checkerboard_decomposition(L):
        """ Recursively loops through the dimensions and returns an mask for the
            sublattices A & B defined as

                |    |    |    |
                B -- A -- B -- A --
                |    |    |    |
                A -- B -- A -- B --

            Convention: sublattice 0 (A) starts at the lower left site.
        """
        if len(L) == 1:
            return [0, 1] * (L[0]//2)
        b = GaussLattice._checkerboard_decomposition(L[:-1])
        return (b + [1-c for c in b]) * (L[-1]//2)


    def checkerboard_indices(self):
        """ Returns the indidces of the sublattices.
        """
        l = GaussLattice._checkerboard_decomposition(self.L)
        return [i for i, s in enumerate(l) if not s], [i for i, s in enumerate(l) if s]


    def base_to_link(self, base_str):
        """ Takes in a list of basis indicies and and returns the link occupation
            representation, i.e., an integer whose binary representation reflects
            occupied/unoccupied links

            Example: A 2x2 system has 8 links such that the maps read

                [5, 5] ---> 11111111
                [0, 0] ---> 00000000

            where the basis vertex '5' is the full vertex and '0' is the empty
            vertex. Intermediate ones work accordingly.

        """
        link_state = 0
        for i, bi in enumerate(base_str):
            link_state += self.lattice_base[i][bi]
        return link_state


    def make_charge_background(self, charges):
        """ Creates all local constraints for the GL, including potentially
            occuring static charges (they change the local constraints on a
            given vertex).
        """
        if len(charges) != 2:
            raise ValueError("Wrong charge configuration specified.")

        cu, cd = charges
        if len(cu) != len(cd):
            raise ValueError("Positive and negative charges need to cancel.")

        # Place the charges on the lattice.
        gl_constraints = [0]*self.S[-1]
        for c in cu:
            gl_constraints[c] = 1
        for c in cd:
            gl_constraints[c] = -1
        return gl_constraints


    def get_vertex_links(self, i):
        """ Returns the indices [+x, -x, +y, -y, ...] of the links for the i-th
            vertex *on the full lattice* (not the sublattice) under consideration
            of periodic boundary conditions.

            Note: works only for L^d lattices for now.
        """
        # Shorthand to avoid self all the time.
        L, d, S = self.L, self.d, self.S

        # Index in bit string.
        j = d*i

        # Loop through the dimensions and add the indicies of + and - directions.
        ind = []
        for k in range(1,d+1):
            # Step forward is always 'on site'.
            ind.append(j+k-1)

            # Sep backward in k direction.
            l = j + k -1 - d*S[k-1]
            if i%S[k] < S[k-1]:
                l = l + d*S[k]
            ind.append(l)
        return ind


    def set_vertex_links(self, i, basis_vertex):
        """ Sets the i-th vertex to be the basis_vertex, which is given as
            occupation numbers in the form

                basic_vertex = [n(+x), n(-x), n(+y), n(-y), ... ]

            Returns a number which represents an empty lattice with only the
            i-th vertex set to the specified basis_vertex.
        """
        d = len(basis_vertex)//2
        lvert = self.get_vertex_links(i)
        b = 0
        for k in range(2*d):
            b += (1<<lvert[k]) * basis_vertex[k]
        return b


    def find_vertex_base(self):
        """ Returns the basis states which reflect the allowed vertices on the
            d-dimensional lattice that obey the Gauss Law.

            Expected number of basis states:
                general:   (2*d!)/(d!*d!)
                2D:        6
                3D:        20

            Not the smartest algorithm (loops through all states) but it does
            the trick.
        """
        d = self.d

        # Produces all possible combinations.
        all_vertices = product([0,1], repeat=2*d)

        # Loops throgh and checks if the GL is valid, if so, adds it to permissible
        # basis states.
        base = {k:[] for k in range(-d,d+1,1)}
        for vertex in all_vertices:
            s = 0
            for i in range(d):
                s += (vertex[2*i] - vertex[2*i+1])
            base[s].append(vertex)
        return base


    def create_spatial_base_elements(self, ind_bs, static_charges, vertex_base):
        """ Creates the basis elements at every specified lattice index in ind_bs.
            This is required because in the presence of static charges the allowed
            vertices are not the same for all sites. This is the most basic version
            which holds an array of allowed vertices for every lattice site - it
            is also the most general way.

            Notes:
                -   This is the only place where the internal representation of
                    the states is defined - change it here and it works differently
                    (the only thing is that those things need to be able to support
                    appending by sum, something like a string or a list).

        """
        spatial_base_elements = []
        for i, bs in enumerate(ind_bs):
            spatial_base_elements.append([[k] for k in range(len(vertex_base[static_charges[bs]]))])
        return spatial_base_elements


    def construct_lattice_basis(self, lattice_vertices, static_charges, vertex_base):
        """ Constructs the basis states (which are encoded in the array base) on every
            (sub)lattice site. This has the benefit that periodic boundary conditions are
            directly implemented and the basis states can directly be looked up in a 2D
            array.

            Ordering: linear, starting bottom left.
        """
        bg = []
        for i, bs in enumerate(lattice_vertices):
            site_basis = []
            for bvert in vertex_base[static_charges[bs]]:
                site_basis.append(self.set_vertex_links(bs, bvert))
            bg.append(site_basis)
        return bg

    def create_gls_masks(self, gls):
        """ Returns masks to check the violation of Gauss law. There are two
            contributions: positive and negative. Both are returned as separate
            dictioinaries that map the lattice index (on the full lattice) to the
            indices and the integer representation of the lattice with only the
            two corresponding indicies set.
        """
        mpos, mneg = {}, {}
        for j, s in enumerate(gls):
            i = self.get_vertex_links(s)
            mpos[s] = i[::2]
            mneg[s] = i[1::2]

            platt = 0
            for p in mpos[s]:
                platt += 1 << p
            mpos[s].append(platt)

            nlatt = 0
            for n in mneg[s]:
                nlatt += 1 << n
            mneg[s].append(nlatt)

        return mpos, mneg


    def find_checkable_gls(self, bs, gls):
        """ Returns a list of lists, each of which hold the indices of the GLS that
            can be checked with partial information of the length of the index. I.e.,
            the i-th list will hold all states that are checkable with a basis string
            of length i.

            Checkable referst to the ability to check whether the lattice configuration
            obeys the Gauss Law from partial information. This is needed for a recursive
            exploration of the state space (which then essentially cuts out entire sub-
            trees and therefore speeds up the search considerably).

            The strategy is to put the full vertex (that is the vertex with full fermion
            occupation at all links) at the first basis site (BS) and check if any Gauss-
            Law-Site (GLS) is fully occupied. If yes, it is added to the checkable sites
            at this order. Then put the full vertex on the second BS and so forth until
            all BS are exhausted.

            Could be optimized, but it does the trick (in arbitrary dimension).
        """
        d = self.d

        # Create the full vertex.
        full_vertex = [1,1]*d
        cstates = []
        latt = 0
        for b in bs:
            # temp = self.set_vertex_links(b, full_vertex)
            latt += self.set_vertex_links(b, full_vertex)
            states = []
            for g in gls:
                # Count if the state is fully surrounded (i.e. if all four adjacent
                # lattice sites are occupied).
                if sum([latt >> i & 1 for i in self.get_vertex_links(g)]) == 2*d:
                    states.append(g)
            cstates.append(states)

        return cstates


    def check_lattice(self, latt, level):
        """ Checks whether the provided lattice string is valid. Also works on
            partial lattice information, however, only check the "checkable"
            states in this case.
        """
        # Loop through all checkable sates at this order.
        for g in self.checkable_gls[level]:
            # We initialize with the static charge, then add the field lines,
            # in the end this should add up to 0 for GL to be valid.
            s = -self.static_charges[g]

            # Sum positive contributions.
            platt = latt & self.mpos[g][-1]
            for p in self.mpos[g][:-1]:
                s += (platt >> p) & 1

            # Sum negative contributions.
            nlatt = latt & self.mneg[g][-1]
            for n in self.mneg[g][:-1]:
                s -= (nlatt >> n) & 1

            # If nonzero, GL is violated.
            if s:
                return False
        return True


    def _find_states(self, prefix):
        """ Finds all possible states of length L that can be constructed from the base,
            which holds all possible basis states.
        """
        # Convert to binary representation.
        latt = self.base_to_link(prefix)
        l = len(prefix)

        # For any system with static charges the condition with the length of the prefix
        # is crucial!
        if not l or self.check_lattice(latt, len(prefix)-1):
            if not l-self.N_sublattice:
                self._collect_state(latt)
                return
            for b in self.spatial_base_elements[l]:
                self._find_states(prefix+b)
            return
        return


    def find_states(self):
        """ Wraps the recursive function for external use.
        """
        self._find_states([])
        if self.write_states:
            self._flush_buffer()
        return self.winding_bins


    def _collect_state(self, latt):
        """ Does all the counting and whatever else is needed.
        """
        # Optional output for 3D calculations, so that we see some progress.
        s = self.winding_bins.sum()
        if not s % 100000:
            print(datetime.datetime.now().strftime("%H:%M:%S") + ' - {:d}'.format(s))

        if self.use_winding:
            w = self.get_winding_numbers(latt)
            self.winding_bins[tuple(w)] += 1

            if self.write_states:
                self.buffer.put([latt] + w)
        else:
            self.winding_bins[0] += 1
            if self.write_states:
                self.buffer.put([latt, self.ds_label])

        # Flush buffer, if full.
        if self.buffer.full():
            self._flush_buffer()


    def get_winding_numbers(self, latt):
        """ Computes the winding numbers for every direction from a state in the
            basis-vertex string representation.
        """
        # Compute all Winding numbers.
        w = []
        for i, wm in enumerate(self.winding_masks):
            winding = 0
            for k in wm:
                winding += (latt >> int(k))&1
            # w[i] = w[i] // self.L[i]
            w.append(winding)
        return w


    def prepare_winding_numbers(self):
        """ Computes the indices to be summed over for the different winding
            numbers.
        """
        L, S, d = self.L, self.S, self.d

        # Find the indices and shifts along the axes. The indices in this case
        # are of the site on the lattice, without pointing to a specific link yet.
        # The shift moves a piont on the axis to the opposite side of the lattice.
        e = []
        for i in range(1,d+1):
            e.append(np.arange(0,S[i],S[i-1])*d)

        if self.d == 2:
            winding_bins = np.zeros(shape=tuple(L[::-1]+1), dtype=np.int)

            # In 2D these are just along one axis.
            winding_masks = [
                e[1],
                e[0]+1
            ]

        if self.d == 3:
            shape = tuple(np.array([
                L[1]*L[2] + 1,
                L[0]*L[2] + 1,
                L[0]*L[1] + 1
            ]))
            winding_bins = np.zeros(shape=shape, dtype=np.int)

            winding_masks = []

            # In 3D we have to sum over entire faces of the cube - here we
            # construct the indices for all of them.
            wmx = []
            for j in range(L[2]):
                wmx = np.concatenate((wmx, np.arange(j*3*S[2], 3*(j*S[2]+S[2]), 3*S[1])))
            winding_masks.append(np.array(wmx, dtype=np.int))
            # winding_masks.append(np.arange(0,3*S[3],3))


            wmy = []
            for j in range(L[2]):
                wmy = np.concatenate((wmy, np.arange(j*3*S[2], 3*(j*S[2]+S[1]), 3*S[0])))
            winding_masks.append(np.array(wmy, dtype=np.int)+1)
            # winding_masks.append(np.arange(1,3*S[3],3))

            wmz = []
            for j in range(L[0]):
                wmz = np.concatenate((wmz, np.arange(j*3*S[1], 3*(j*S[1]+S[1]), 3*S[0])))
            winding_masks.append(np.array(wmz, dtype=np.int)+2)
            # winding_masks.append(np.arange(2,3*S[3],3))

        return winding_bins, winding_masks



    # ==========================================================================
    # I/O stuff.

    def _init_file(self, state_file, append=False):
        """ Sets up the state output.
        """
        self.state_file = self.out_dir + '/' + state_file
        self.output_format = state_file.split('.')[-1]

        if self.output_format == 'hdf5':
            with hdf.File(self.state_file, 'a' if append else 'w') as f:
                # We can loop through all winding number sectors with the product
                # functions, which is essentially a cartesian product generator.
                ds_labels = []
                if self.use_winding:
                    winding_indicies = product(*map(lambda n: range(n), self.winding_bins.shape))
                    for ws in winding_indicies:
                        ds_labels.append(winding_tag(ws))
                else:
                    ds_labels = [self.ds_label]

                for ds_label in ds_labels:
                    if ds_label in f:
                        del f[ds_label]
                    dset = f.create_dataset(
                        ds_label,
                        (0,),
                        maxshape=(None,),
                        dtype='i8',
                        chunks=True
                    )
        else:
            # Attention: truncates existing file.
            with open(self.state_file, 'w') as f:
                # This is a dimensional "limitation" - works only for up to 3D.
                labels = np.array(['x', 'y', 'z'])[:self.d]
                f.write('state,' + ('w_{:s},'*self.d).format(*labels)[:-1] + '\n')


    def _flush_buffer(self):
        """ Writes the buffer to file and clears it.
        """
        # Writing with the HDF5 format requires some extra work, since we want to
        # also sort it by winding sector. To do this, we need to first sort the
        # corresponding winding sectors and then loop through all of them to append
        # to the file.
        # We could, alternatively, append to the HDF5 datasets one-by-one, but
        # this will likely increase the overhead in file storage time and disk
        # space significantly (I say that without testing it though).
        if self.output_format == 'hdf5':

            # Sorting (probably inefficient, but this is a rare task).
            wn_dict = {}
            try:
                for line in iter(self.buffer.get_nowait, None):
                    state, ws = line[0], tuple(line[1:])
                    wlist = wn_dict.get(ws)
                    if not wlist:
                        wn_dict[ws] = [state]
                    else:
                        # For all the C-programmers who may or may not make it to this
                        # point: this works, because Python defaults to call-by-reference!
                        wlist.append(state)
            except Empty:
                pass

            # Write to HDF5 file.
            with hdf.File(self.state_file, 'a') as f:
                for ws, states in wn_dict.items():
                    ds_tag = winding_tag(ws) if self.use_winding else ws[0] # Somewhat dirty hack to make the code work for non-winding sorted stuff as well.
                    dset = f[ds_tag]
                    dset.resize(dset.shape[0]+len(states), axis=0)
                    dset[-len(states):] = states

        else:
            # Iterates through the queue and empties it in a FIFO manner.
            with open(self.state_file, 'a+') as f:
                try:
                    for line in iter(self.buffer.get_nowait, None):
                        f.write(','.join(map(str, line))+'\n')
                except Empty:
                    pass
