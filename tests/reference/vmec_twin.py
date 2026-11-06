"""Run VMEC++ on a VMEC INDATA file (PB3D's --export_HEL output) and save the wout file.

Usage: python vmec_twin.py <indata-file> <wout-file>   (python with the vmecpp package)
The resolution is set here: five multigrid stages up to 319 surfaces, force
tolerance 1e-15. Prints q on axis and at the edge, the boundary coefficients and
the force residuals, for the comparison with the HELENA file the input came from.
"""
import sys
import numpy as np
import vmecpp
import netCDF4

inp, wout = sys.argv[1], sys.argv[2]
vi = vmecpp.VmecInput.from_file(inp)
upd = {"ns_array": np.array([19, 39, 79, 159, 319]),
       "ftol_array": np.array([1e-8, 1e-10, 1e-12, 1e-14, 1e-15]),
       "niter_array": np.array([20000] * 5),
       "return_outputs_even_if_not_converged": True}
vi = vi.model_copy(update={k: v for k, v in upd.items() if hasattr(vi, k)})
print("mpol", vi.mpol, "ntor", vi.ntor, "ncurr", vi.ncurr, "phiedge", vi.phiedge, "piota_type", vi.piota_type)
out = vmecpp.run(vi, max_threads=8, verbose=False)
out.wout.save(wout)

d = netCDF4.Dataset(wout)
g = lambda k: d.variables[k][:]
iot = g("iotaf")
print("ns", int(g("ns")), "mpol", int(g("mpol")),
      "iota axis %.6f edge %.6f -> q axis %.5f edge %.5f" % (iot[0], iot[-1], 1 / iot[0], 1 / iot[-1]))
print("phi edge", float(g("phi")[-1]), "Raxis", float(g("raxis_cc")[0]), "aspect", float(g("aspect")),
      "betatotal", float(g("betatotal")), "volume", float(g("volume_p")))
rmnc = g("rmnc"); zmns = g("zmns")
print("edge R(m=0..3)", [float(rmnc[-1, i]) for i in range(4)], "edge Z(m=1)", float(zmns[-1, 1]))
for k in ("fsqr", "fsqz", "fsql", "ier_flag", "ftolv"):
    if k in d.variables:
        print(k, float(g(k)))
