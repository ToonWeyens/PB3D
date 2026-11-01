"""Write a HELENA input for the cbm18a plasma cut at the pedestal foot.
Usage: make_cut_input.py <extended fort.10> <out fort.10> B XIAB coarse|fine [NPTS_cut=217]"""
import sys, re, numpy as np
src, dst, bpar, xiab, res = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), sys.argv[5]
ncut = int(sys.argv[6]) if len(sys.argv) > 6 else 217
txt = open(src).read()
dpr = re.findall(r'dpr\( *\d+\)= *([-+0-9.eE]+)', txt); df2 = re.findall(r'df2\( *\d+\)= *([-+0-9.eE]+)', txt)
am = np.load('foot_am.npy'); tg = np.load('foot_targets.npy'); eps = tg[2]
mfm = 24
L = [' &SHAPE IAS = 0, ISHAPE = 2, MFM = %d, MHARM = 128, ELLIP=1.0, TRIA=0.0, QUAD=0.,' % mfm]
L += ['   FM(%d) = %.10e,' % (m + 1, am[m]) for m in range(mfm)]
L += [' &END', ' &PROFILE', '    IPAI=11, EPI=1.0, FPI=1.0, GPI=0.01', '    IGAM=11, EGA=1.0, FGA=1.0, GGA = 0.01,', '    ICUR=0, ECUR=1.0, FCUR=1.0, GCUR=0.01,', '    NPTS = %d' % ncut]
L += [' dpr(%3d)=%s, df2(%3d)=%s,' % (i + 1, dpr[i], i + 1, df2[i]) for i in range(ncut)]
L += [' &END', ' &PHYS', '  EPS   = %.8f' % eps, '  ALFA  = 5.5', '  B     = %.10e' % bpar, '  BETAP = -1.0', '  XIAB  = %.6f' % xiab, '  RVAC  = 3.0', '  BVAC  = 2.0', ' &END', ' &NUM']
if res == 'coarse':
    L += ['    NR = 101, NP = 65, NRMAP = 201, NPMAP = 129, NCHI = 129', '    NRCUR = 101, NPCUR = 65, ERRCUR = 1.e-6, NITER = 100, NMESH = 100']
else:
    L += ['    NR = 251, NP = 201, NRMAP = 401, NPMAP = 513, NCHI = 513', '    NRCUR = 251, NPCUR = 128, ERRCUR = 1.e-6, NITER = 100, NMESH = 100']
L += [' &END', ' &PRI  NPR1=1 &END', ' &PLOT NPL1=1 &END', ' &BALL        &END', '']
open(dst, 'w').write('\n'.join(L))
