""" ----------------------------------------------------------------------------

    lattice_object.py - LR, July 2020

    This defines a class that represents a lattice and allows to perform several
    symmetry operations on it. This is explicitly not for computation, as it is
    a quite inefficient representation. But this can be used to investigate
    the effect of certain operators or rotation and to figure out the symmetries
    of the underlying problem.

    Several old snippets of code have been mashed together to achive this. The
    code could be much cleaner and more structured, however, it should ge the
    stuff done.

---------------------------------------------------------------------------- """
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as a3
import matplotlib.colors as colors
import matplotlib.patches as mpatches
import mpl_toolkits.mplot3d.art3d as art3d
import matplotlib.patheffects as patheffects
import matplotlib.patches as patches


# We need a builder object so that we can perform some operations more easily.
import sys
sys.path.append('../python_gauss_lattice/')
from gauss_lattice import HamiltonianBuilder


class Vertex(object):
    """ Represents a vertex in the grid with all its neighbours (previous and next).
    """
    def __init__(self, links):
        """ Takes in the links in the form [x, -x, y, -y, z, -z] and stores them.
        """
        self.v = np.zeros(shape=(3,3,3), dtype=int)
        self.v[2,1,1] = links[0]
        self.v[0,1,1] = links[1]

        self.v[1,2,1] = links[2]
        self.v[1,0,1] = links[3]

        #self.v[1,1,2] = links[4]
        #self.v[1,1,0] = links[5]

    def bit_string(self):
        """ Returns the partial bit string, i.e., integer for the vertex.
            Only the "upper" entries are of importance here.
        """
        return self.v[2,1,1] + 2*self.v[1,2,1] + 4*self.v[1,1,2]

    def flip(self):
        self.v[2,1,1] = 1 - self.v[2,1,1]
        self.v[0,1,1] = 1 - self.v[0,1,1]
        self.v[1,2,1] = 1 - self.v[1,2,1]
        self.v[1,0,1] = 1 - self.v[1,0,1]
        self.v[1,1,2] = 1 - self.v[1,1,2]
        self.v[1,1,0] = 1 - self.v[1,1,0]

    def rot90(self, axes, k=1):
        self.v = np.rot90(self.v, k=k, axes=axes)

    def __repr__(self):
        return "{:03b}".format(self.bit_string())


class LatticeObject2D(object):
    """ Represents a lattice and lets us perform some operations on it.

        Quite inefficient when thinking of acutal memory consumption, but this is just for
        the purpose of symmetry exploration - way too slow for actual computations.
    """
    def __init__(self, state, L, quiet=False):
        """ Takes an integer as the bit representation for the lattice.
        """
        if not quiet:
            print('Setting up lattice {:d}'.format(state))
        self._update(L)

        self.vertices = np.zeros(shape=L, dtype=Vertex)
        for i in range(self.S[-1]):
            link_indicies = self.get_vertex_links(i)
            links = [(state>>k)&1 for k in link_indicies]
            x, y = self.get_2d_indicies(i)
            self.vertices[x,y] = Vertex(links)

    def _update(self, L):
        self.L = L
        self.d = len(L)
        self.S = [1]
        for l in L:
            self.S.append(l*self.S[-1])
        self.nb = self.S[-1]*self.d


    def get_3d_indicies(self, i):
        """ Takes in a linear index and maps it to the 3D array.
        """
        return i%self.S[1], i%self.S[2]//self.S[1], i%self.S[3]//self.S[2]

    def get_2d_indicies(self, i):
        """ Takes in a linear index and maps it to the 3D array.
        """
        return i%self.S[1], i%self.S[2]//self.S[1]

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

    def to_int(self):
        """ Get back the integer representation.

            Note: Here it is better to explicitly loop through the indicies to keep
            the correct ordering.
        """
        state = 0
        for i in range(self.S[-1]):
            x, y = self.get_2d_indicies(i)
            state += 2**(self.d*i)*int(self.vertices[x,y].bit_string())
        return int(state)

    def to_bin(self):
        return ('|{:0'+str(self.nb)+'b}> ({:d})').format(self.to_int(), self.to_int())

    # ----
    # Drawing stuff.

    def _get_drawing_links(self, b, a):
        """ Takes coordinates and returns the links to a vertex.
        """
        return [(*b[:k], b[k]+a,*b[k+1:]) for k in range(len(b))]


    def _particle_coordinates(self, link, a=2):
        """ Maps a link to the coordinates on the plot. Returns start and endpoint.
        """
        s = link // len(self.L)
        d = link % len(self.L)

        base = [
            a * (s %self. L[0]),
            a * ((s // self.L[0]) % self.L[1]),
            a * ((s // (self.L[0]*self.L[1])) % self.L[2])
        ]
        l = list(self._get_drawing_links(base, a)[d])
        l[d] -= a/2
        return l


    def _link_coordinates(self, link, a=2):
        """ Maps a link to the coordinates on the plot. Returns start and endpoint.
        """
        s = link // len(self.L)
        d = link % len(self.L)

        base = [
            a * (s % self.L[0]),
            a * ((s // self.L[0]) % self.L[1])
        ]
        l = self._get_drawing_links(base, a)
        return base, l[d]


    def _draw_2D_lattice(self, a=2.0, pbc=True):
        all_links = []
        for y in range(self.L[1]):
            for x in range(self.L[0]):
                links = self._get_drawing_links([x*a, y*a], a)
                for l in links:
                    all_links.append([[x*a, l[0]], [y*a, l[1]]])
                    if (len(all_links)-1) == 100:
                        self.ax.plot(*all_links[-1], color='red')
                    else:
                        self.ax.plot(*all_links[-1],
                                    color='black',
                                    # color='#a3a2a2',
                                    lw=1.5,
                                    zorder=10-a*y
                            )

        return all_links

    def _draw_particle(self, pos, occupied, a=2.0, colors=['#416780', '#5C60C0'], type='particle'):
        """ Draws a particle.
        """
        # Direction of link of particle.
        dir = pos%len(self.dlinks[0])

        # Particle coordinate from link coordinates (left/lower corner).
        coord = list(map(lambda x: [x[0]], self.dlinks[pos]))
        coord[dir][0] += a/2 # Shift in appropriate direction.

        if (type == 'particle') and occupied:
            self.ax.plot(coord[0], coord[1],
                    marker='o',
                    markeredgecolor=colors[0],
                    markerfacecolor=colors[1],
                    markersize=8,
                    markeredgewidth=1,
                    zorder=100
            )
        elif type =='spin':

            self.ax.plot(coord[0], coord[1],
                    marker=('^' if dir else '>') if occupied else ('v' if dir else '<'),
                    markeredgecolor=colors[0],
                    markerfacecolor=colors[1],
                    markersize=8,
                    markeredgewidth=1,
                    zorder=100
            )

    def _draw_state(self, colors=['black', '#5C60C0'], type='particle'):
        state = self.to_int()
        for k in range(self.nb):
            self._draw_particle(k, (state>>k)&1, colors=colors, type=type)

    def _highlight_plaquette(self, p, color='gray', alpha=0.2, a=2):
        lcs = [self._link_coordinates(v) for v in p]
        base = lcs[0][0]
        rect = patches.Rectangle(
            base, a, a,
            linewidth=1, edgecolor='none', facecolor=color, alpha=alpha
        )
        self.ax.add_patch(rect)

    def _find_flippable_plaquettes(self, builder):
        """ Finds the flippable plaquettes and returns them.
        """
        state = self.to_int()
        pflips = []
        for p in builder.plaquettes:
            if builder.apply_u_dagger(state, p)[0] or builder.apply_u(state, p)[0]:
                pflips.append(p)
        return pflips


    def draw(self, a=2, axis=None, label=None, plaquettes=True, type='particle'):
        """ 3D plot for the current lattice.
        """
        if axis is None:
            self.fig = plt.figure()
            self.fig.set_size_inches(2.6,2.2)
            self.ax = self.fig.add_subplot(111)
        else:
            self.ax = axis


        # Draw the actual system.
        self.dlinks = self._draw_2D_lattice()
        points = self._draw_state(type=type)

        # Get flippable plaquettes & color them.
        if plaquettes:
            builder = HamiltonianBuilder({'L':self.L}, [], silent=True)
            pflip = self._find_flippable_plaquettes(builder)
            print(f'# of flippable plaquettes: {len(pflip)}')
            for p in builder.plaquettes:
                if p in pflip:
                    pass
                    self._highlight_plaquette(p[:-1], color='green', alpha=0.1)
                else:
                    pass
                    self._highlight_plaquette(p[:-1], color='red', alpha=0.1)

        # Set viewpoint correctly.
        self.ax.set_axis_off()
        # self.ax.set_xlim(0, a*self.L[0])
        # self.ax.set_ylim(0, a*self.L[1])

        if label is not None:
            self.ax.text(-0.5, 0, self.ax.get_zlim()[1]+0.5, label, fontsize=24)


    # -----
    # Operators.

    def apply_charge_conjugation(self):
        """ Flips all occupancies (spins) on the lattice.
        """
        for f in self.vertices.flat:
            f.flip()

    def apply_parity_flip(self):
        """ Reverses all axes.
        """
        for k in range(self.d):
            self.vertices = np.flip(self.vertices, axis=k)

    def apply_translation(self, axis, extent=1):
        """ Translates steps along a specified axis.
        """
        self.vertices = np.roll(self.vertices, shift=extent, axis=axis)
