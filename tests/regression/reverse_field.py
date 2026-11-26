#!/usr/bin/env python3
"""Write a copy of a VMEC wout file with the magnetic field reversed (issue #19).

With B replaced by -B, and with it the current J by -J, the pressure balance
and the flux surfaces stay the same, and so does the ideal MHD stability. In
the wout file, the geometry (R, Z, lambda, the Jacobian), the rotational
transform, |B| and the pressure stay the same; the fluxes and the components
of B and J change sign. The list below holds the variables of that kind that
a wout file can have (VMEC++ and VMEC2000); the names that end in the other
letter (c for s, s for c) exist only in files without stellarator symmetry.
VMEC++ solving the same input with PHIEDGE and CURTOR of the opposite sign
gives exactly this copy.

Usage: reverse_field.py wout_in.nc wout_out.nc
"""
import shutil
import sys

import numpy as np

OPPOSITE = ['phi', 'phipf', 'phips', 'chi', 'chipf',
            'bsubsmns', 'bsubumnc', 'bsubvmnc', 'bsupumnc', 'bsupvmnc',
            'bsubsmnc', 'bsubumns', 'bsubvmns', 'bsupumns', 'bsupvmns',
            'buco', 'bvco', 'ctor', 'currumnc', 'currvmnc', 'currumns', 'currvmns',
            'jcuru', 'jcurv', 'rbtor', 'rbtor0', 'b0', 'bdotgradv']


def reverse(path):
    """Change the sign of the variables of OPPOSITE that the file has, in place."""
    done = []
    try:
        import netCDF4
        with netCDF4.Dataset(path, 'r+') as d:
            for k in OPPOSITE:
                if k in d.variables:
                    d.variables[k][...] = -d.variables[k][...]
                    done.append(k)
        return done
    except ImportError:
        pass
    import h5py
    if h5py.is_hdf5(path):                                                      # netCDF 4 (e.g. VMEC++)
        with h5py.File(path, 'r+') as d:
            for k in OPPOSITE:
                if k in d:
                    d[k][...] = -d[k][...]
                    done.append(k)
    else:                                                                       # netCDF 3 (e.g. VMEC2000)
        from scipy.io import netcdf_file
        with netcdf_file(path, 'a', mmap=False) as d:
            for k in OPPOSITE:
                if k in d.variables:
                    d.variables[k][...] = -np.array(d.variables[k][...])
                    done.append(k)
    return done


def main():
    src, dst = sys.argv[1], sys.argv[2]
    shutil.copyfile(src, dst)
    done = reverse(dst)
    print('sign changed in %d variables: %s' % (len(done), ' '.join(done)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
