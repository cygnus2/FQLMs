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

        self.v[1,1,2] = links[4]
        self.v[1,1,0] = links[5]

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


class LatticeObject(object):
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
            x, y, z = self.get_3d_indicies(i)
            self.vertices[x,y,z] = Vertex(links)

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
            x, y, z = self.get_3d_indicies(i)
            state += 2**(3*i)*self.vertices[x,y,z].bit_string()
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
            a * ((s // self.L[0]) % self.L[1]),
            a * ((s // (self.L[0]*self.L[1])) % self.L[2])
        ]
        l = self._get_drawing_links(base, a)
        return base, l[d]


    def _draw_3D_lattice(self, a=2.0, pbc=True):
        all_links = []
        for z in range(self.L[2]):
            for y in range(self.L[1]):
                for x in range(self.L[0]):
                    links = self._get_drawing_links([x*a, y*a, z*a], a)
                    for l in links:
                        all_links.append([[x*a, l[0]], [y*a, l[1]], [z*a, l[2]]])
                        if (len(all_links)-1) == 100:
                            self.ax.plot(*all_links[-1], color='red')
                        else:
                            self.ax.plot(*all_links[-1],
                                    color='#838282',
                                    lw=2.5,
                                    zorder=10-a*y
                            )

        return all_links

    def _draw_particle(self, pos, links, a=2.0, colors=['#FCC148', '#F9D792']):
        coord = list(map(lambda x: [x[0]], links[pos]))
        coord[pos%len(links[0])][0] += a/2
        self.ax.plot(xs=coord[0], ys=coord[1], zs=coord[2],
                marker='o',
                markeredgecolor=colors[0],
                markerfacecolor=colors[1],
                markersize=10,
                markeredgewidth=2,
                zorder=10-a*(coord[1][0]//a)+1
        )


    def _draw_state(self, links, colors=['#FCC148', '#F9D792']):
        state = self.to_int()
        for k in range(self.nb):
            if (state>>k)&1:
                self._draw_particle(k, links, colors=colors)

    def _highlight_plaquette(self, p, color='gray', alpha=0.2, a=2):
        lcs = [self._link_coordinates(v) for v in p]

        base = lcs[0][0]
        shift_1 = np.array(lcs[0][1]) - np.array(lcs[0][0])
        shift_2 = np.array(lcs[1][1]) - np.array(lcs[1][0])
        shift_3 = np.array(lcs[2][1]) - np.array(lcs[2][0])
        pts = np.array([
            base,
            base + shift_1,
            base + shift_1 + shift_2,
            base - shift_1 + shift_2 + shift_3
        ])
        tri = a3.art3d.Poly3DCollection([pts], alpha=alpha, color=color)
        self.ax.add_collection3d(tri)

    def _find_flippable_plaquettes(self, builder):
        """ Finds the flippable plaquettes and returns them.
        """
        state = self.to_int()
        pflips = []
        for p in builder.plaquettes:
            if builder.apply_u_dagger(state, p)[0] or builder.apply_u(state, p)[0]:
                pflips.append(p)
        return pflips

    def draw(self, a=2):
        """ 3D plot for the current lattice.
        """
        with plt.style.context('seaborn-notebook'):
            self.fig = plt.figure()
            self.fig.set_size_inches(8,6)
            self.ax = self.fig.add_subplot(111, projection='3d')

            self.ax.view_init(elev=20, azim=-75)
            self.ax.set_axis_off()

            # Draw the actual system.
            links = self._draw_3D_lattice()
            points = self._draw_state(links)

            # Get flippable plaquettes & color them.
            builder = HamiltonianBuilder({'L':self.L}, [], silent=True)
            pflip = self._find_flippable_plaquettes(builder)
            print(f'# of flippable plaquettes: {len(pflip)}')
            for p in builder.plaquettes:
                if p in pflip:
                    pass
#                     self._highlight_plaquette(p[:-1], color='green', alpha=0.1)
                else:
                    pass
                    self._highlight_plaquette(p[:-1], color='red', alpha=0.1)

                    
                    
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
        
    def _rotation_sequence(self, rotations, steps):
        """ Performs a sequence of rotations.
        """
        for i in range(steps):
            for ax, steps in rotations:
                self.apply_c4(ax, steps=steps)
                
    def apply_c2(self, axis, steps=1):
        """ Rotation about axes on oposite edges.
            
            Notes:
             - 6 axes with 1 rotation of 180 degrees = 1 step
             - 2 steps are the identity
             - 6 elements
             
             Axis numeration, denotes connected links:
              0 : (2,4) - (5,7)
              1 : (1,3) - (6,8)
              2 : (5,6) - (3,4)
              3 : (1,2) - (7,8)
              4 : (2,6) - (3,7)
              5 : (1,5) - (4,8)
        """
        rotations = [
            [['z', 2], ['y', 1]],
            [['z', 2], ['y', -1]],
            [['y', 2], ['x', 1]],
            [['y', 2], ['x', -1]],
            [['x', 2], ['z', 1]],
            [['x', 2], ['z', -1]],
        ]
        self._rotation_sequence(rotations[axis], steps)
    
    def apply_c4(self, axis, steps=1):
        """ Rotation around coordinate axes.  
            
            Notes:
             - 3 axes with 3 rotations of 90 degrees = 1 step
             - 4 steps are the identity
             - 9 elements
        """
        ax = {
            'x' : (1,2),
            'y' : (2,0),
            'z' : (0,1), 
        }
        self.vertices = np.rot90(self.vertices, k=steps, axes=ax[axis])
        for f in self.vertices.flat:
            f.rot90(k=steps, axes=ax[axis])
        
        # Needed for unequal lattices - to get the indicies right.
        self._update(self.vertices.shape)
    
    def apply_c3(self, axis, steps=1):
        """ Rotation around the diagonals between to opposite corners. 
            
            Notes:
             - 4 axes with 2 rotations of 120 degrees = 1 step
             - 3 steps are the identity
             - 8 elements
             
             Axis numeration, denotes connected corners:
              0 : (1) - (8)
              1 : (2) - (7)
              2 : (4) - (5)
              3 : (3) - (6)
        """
        rotations = [
            [['y', -1], ['z', -1]],
            [['y', 1], ['z', 1]],
            [['y', 1], ['z', -1]],
            [['z', -1], ['y', 1]],
        ]
        self._rotation_sequence(rotations[axis], steps)