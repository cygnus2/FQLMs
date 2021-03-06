""" ----------------------------------------------------------------------------

    ed_result_wrapper.py - LR, January 2021

    A class that represents results for a given HDF5 output of the code. This
    provides some convenience when plotting stuff.

    Notes:
     -  Works for both, Python and Julia output (since they should be the same)
     -  Naming conventions hard-coded, best if left untouched.

---------------------------------------------------------------------------- """
import h5py as hdf
import numpy as np
import pandas as pd
from copy import copy


class EDResult(object):
    def __init__(self, param, datadir, quiet=False, datafile=None):
        """ Construct the filename and read contents to a DataFrame. Eigenstates
            are only read on demand due to possibly large memory requirements.
        """
        if 'winding_sector' in param:
            self.ws_label = EDResult._winding_tag(param['winding_sector'], shift=param['L'])
        else:
            self.ws_label = 'all-ws'

        if 'static_charges' in param:
            self.charge_label = EDResult._charge_tag(param['static_charges'])
        else:
            self.charge_label = None

        # General fields.
        self.defaults = copy(param)
        for k in self.defaults:
            if isinstance(self.defaults[k], list):
                self.defaults[k] = tuple(self.defaults[k])

        ignored = ['lambda']
        for k in ignored:
            if k in self.defaults:
                del self.defaults[k]

        # Toggle low-energy mode.
        self.le = param.get("low_energy_run")

        # These are the unique keys for the datasets.
        self.keys = ["lambda", "flips"] if self.le else ["lambda"]
        self.keys += list(self.defaults.keys())

        # Datafile.
        if not datafile:
            self.datafile = (datadir+"/{:s}results_{:s}_{:s}_" + ("{:d}x"*len(param["L"]))[:-1] + ".hdf5").format(
                "le_" if self.le else "",
                param['gauge_particles'],
                self.ws_label if self.charge_label is None else self.charge_label,
                *param["L"]
            )
        else:
            self.datafile = datadir+'/'+datafile


        self._read_HDF5()
        if not quiet:
            print(self)

    def _read_HDF5(self):
        """ Converts the provided HDF5 file to a useful format.
        """
        # First collect every value in a separate row.
        rows, all_obs = [], set()
        with hdf.File(self.datafile, 'r') as f:
            if self.le:
                for grp in f:
                    defaults = copy(self.defaults)
                    defaults.update({'flips':int(grp[2:])})
                    r, ao = EDResult._convert_group(f[grp], defaults=defaults)
                    rows += r
                    all_obs = all_obs.union(ao)
            else:
                rows, all_obs = EDResult._convert_group(f, self.defaults)

        # With these rows, we create
        df = pd.DataFrame(rows)

        # To have a useful format, we need to combine lines that are for the same
        # value of lambda and N.
        self.ev = df.groupby(self.keys + ["N"], as_index=False).aggregate({o : "first" for o in all_obs})


    @staticmethod
    def _convert_group(grp, defaults={}):
        """ Takes an open HDF5 group, which represents a container of results,
            and returns all eigenvalues in row format.
        """
        rows, all_obs = [], set()
        for ds in grp:
            obs, *_, lam = ds.split('_')
            if obs != "eigenstates":
                all_obs.add(obs)
                for n, val in enumerate(grp[ds]):
                    row = copy(defaults)
                    row.update({
                        'lambda' : float(lam), # Lambda.
                        obs : val, # Actual value.
                        'N' : n # Number, as in N-th eigenvalue from below.
                    })
                    rows.append(row)
        return rows, all_obs

    def _lam_format(self, lam):
        return "_lam_{:.6f}".format(lam)

    def __repr__(self):
        return "ED results for datafile '" + str(self.datafile) + "'"

    def get_eigenstates(self, lam):
        """ Retrieves the eigenstates for a given lambda.
        """
        with hdf.File(self.datafile, 'r') as f:
            es = f['eigenstates'+self._lam_format(lam)][...]
        return es

    def get_eigenvalues(self, cutoff=20):
        return self.ev

    def get_gs(self):
        """ Returns the ground-state energy as function of lambda.
        """
        gsdf = self.ev[self.ev['N']==0][self.keys+['spectrum']]
        gsdf.columns = self.keys + ['gs_energy']
        return gsdf

    def compute_gap(self, gs=None, skip=['winding_sector']):
        """ Adds the gap to the dataframe.
        """
        if gs is None:
            gs = self.get_gs()

        # Sometimes we have to skip a few keys in order to use this functionality.
        # Example: ground state for different winding sectors.
        # Use with care!
        keys = [k for k in self.keys if k not in skip]

        self.ev = pd.merge(self.ev, gs, on=keys)
        self.ev['delta_e'] = self.ev['spectrum'] - self.ev['gs_energy']


    @staticmethod
    def _winding_tag(ws, labels=['x', 'y', 'z'], shift=None):
        """ Returns the naming convention of the winding datasets.

            The shift is a lattice configuration L = [Lx,Ly,Lz], such
            that the HDF5 datasets may be resolved.
        """
        if shift is not None:
            ws = EDResult._winding_shift(shift, ws)
        wtag = ''
        for k in range(len(ws)):
            wtag += 'w{:s}_{:d}-'.format(labels[k], ws[k])
        return wtag[:-1]


    @staticmethod
    def _winding_shift(L, ws):
        """ Maps between the representation of winding numbers.
        """
        if len(L) == 2:
            shift = np.array(L[::-1]) // 2
        elif len(L) == 3:
            shift = np.array([
                L[1]*L[2] // 2,
                L[0]*L[2] // 2,
                L[0]*L[1] // 2
            ])
        else:
            raise NotImplementedError('Dimension not implemented!')
        return ws + shift


    @staticmethod
    def _charge_tag(static_charge_indicies):
        return (
            'p-'+"-".join(map(lambda x: str(x), static_charge_indicies[0])) +
            "_" +
            'n-'+"-".join(map(lambda x: str(x), static_charge_indicies[1]))
        )
