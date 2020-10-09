""" ----------------------------------------------------------------------------

    combine.py - LR, October 2020

    Combines two sets of lattices to find the correct number of GL lattices
    for a larger geometry.

---------------------------------------------------------------------------- """
import numpy as np
import sys
from copy import copy, deepcopy
from itertools import product
from tqdm import tqdm as tbar



sys.path.append("../")
from gauss_lattice.aux import print_2D_state, bin_str, timeit, read_all_states
from gauss_lattice import GaussLattice


class SimpleLattice(object):
    """ Represents a simple lattice (essentially an integer) and gives also
        a multiplication, so that the combination is easier. For simplicity,
        we assume stacking along the last axis.

        Mainly for diagonistic purposes (way too slow otherwise).
    """
    def __init__(self, state, L):
        """ L can either be an integer or the lattice dimensions.
        """
        self.state = state
        if isinstance(L, int):
            self.L = L
        else:
            self.L = int(np.prod(L))*len(L)

    def __mul__(self, other):
        """ Note that the multiplication is not commutative here, since we're
            chaining together two pieces. Here we assume left multiplication.
        """
        # print(self.L)
        new_state =  self.state + (other.state << self.L)
        return SimpleLattice(new_state, self.L+other.L)


    def __repr__(self):
        # return ("({:3d}) [{:0"+str(self.L)+"b}]").format(self.L, self.state)
        return "({:3d}) [".format(self.L) + bin_str(self.state, self.L) + "]"



def check_lattice(latt, indicies, mpos, mneg):
    """ Checks whether the provided lattice is a valid GL lattice.
    """
    # Loop through all checkable sates at this order.
    for g in indicies:
        s = 0

        # Sum positive contributions.
        platt = latt & mpos[g][-1]
        for p in mpos[g][:-1]:
            s += (platt >> p) & 1

        # Sum negative contributions.
        nlatt = latt & mneg[g][-1]
        for n in mneg[g][:-1]:
            s -= (nlatt >> n) & 1

        # If nonzero, GL is violated.
        if s:
            return False
    return True


def flip_bits(state, bits):
    new_state = state
    for b in bits:
        new_state = new_state^(1<<b)
    return new_state

@timeit(logger=None)
def find_lattices_naive(L):
    gl = GaussLattice(L)
    mpos, mneg = gl.create_gls_masks(list(range(np.prod(L))))

    valid_lattices = []
    for i in range(2**(np.prod(L)*len(L))):
        latt = SimpleLattice(i, L)
        if check_lattice(latt.state, latt.L//len(L), mpos, mneg):
            valid_lattices.append(latt.state)

    print(("Found {:d} valid lattices for "+"x".join(map(str, L))+".").format(len(valid_lattices)))
    return valid_lattices


@timeit(logger=None)
def combine_lattices(lattices_a, lattices_b, La, Lb):
    """ Combines the lattices from two lists to get correct GL states for the third.

        Note: assumes that we combine along the last axis.
    """
    # First, procude the information for the larger lattice.
    Lc = list(La)[:-1] + [La[-1]+Lb[-1]]
    gl = GaussLattice(Lc)

    # Not all indicies have to be checked fot GL validity - here we find which
    # ones we do need to check and prepare the masks for them.
    if len(Lc) == 2:
        check_indicies = [0, 1, np.prod(La), np.prod(La)+1]
        flip_indicies = [0, 2]

    elif len(Lc) == 3:
        s = np.prod(La)
        check_indicies = [0, 1, 2, 3, s, s+1, s+2, s+3]
        flip_indicies = [0, 1, 3, 4, 6, 7, 9, 10]
    else:
        raise ValueError("Wrong dimensions specified.")
    mpos, mneg = gl.create_gls_masks(check_indicies)

    # Produce the masks for flipping.
    flips = []
    for i in product(*[[None, i] for i in flip_indicies]):
        l = list(filter(lambda x: x is not None, list(i)))
        s = sum([1<<k for k in l])
        flips.append(s)



    # This shifts the integers such that they can be combined along the last axis.
    shift = int(np.prod(La)*len(La))

    # found_lattices = []
    # for la in lattices_a:
    #     for bits_a in flips:
    #         latta = flip_bits(la, bits_a)
    #         for lb in lattices_b:
    #             for bits_b in flips:
    #                 lattb = flip_bits(lb, bits_b)
    #
    #                 lattc = latta + (lattb << shift)
    #                 if check_lattice(lattc, check_indicies, mpos, mneg):
    #                     found_lattices.append((lattc))

    found_lattices = []
    clist = list(product(range(len(lattices_a)), range(len(lattices_b))))
    print(len(clist))
    c = 0
    # for i, j in tbar(clist):
    for i in range(len(lattices_a)):
        for j in range(len(lattices_b)):
            c += 1
            for fa in flips:
                latta = lattices_a[i]^fa
                for fb in flips:
                    lattb = lattices_b[j]^fb
                    lattc = latta + (lattb << shift)
                    if check_lattice(lattc, check_indicies, mpos, mneg):
                        found_lattices.append((lattc))
            if c%10 == 0:
                print("hey", c)

    found_lattices = set(found_lattices)
    return set(found_lattices), Lc

La = [2,2,2]
# lattices_a = find_lattices_naive(La)
lattices_a = read_all_states(La, basedir='../output/')
print(("Read {:d} states for "+"x".join(map(str, La))+".").format(len(lattices_a)))

Lb = [2,2,2]
if La == Lb:
    lattices_b = copy(lattices_a)
else:
    # lattices_b = find_lattices_naive(Lb)
    lattices_b = read_all_states(Lb, basedir='../output/')
    print(("Read {:d} states for "+"x".join(map(str, Lb))+".").format(len(lattices_b)))

found_lattices, L_new = combine_lattices(lattices_a, lattices_b, La, Lb)
print(("Found {:d} valid lattices for "+"x".join(map(str, L_new))+" by combination.").format(len(found_lattices)))
