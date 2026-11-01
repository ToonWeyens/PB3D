"""Write a TERPSICHORE control-card file (fixed columns, as tprgl0_block_ap.f reads it).
Usage: python terps_data.py out.data title MM NMIN NMAX MMS NSMIN NSMAX NJ NK IVAC NOWALL AWALL MODELK NSTA AL0 [m_lo m_hi [MPINIT [LPRESS [LCURRF]]]]
Boozer table: all ones for m <= MM, every N; stability table: ones for m_lo..m_hi, every N in NSMIN..NSMAX.
"""
import sys
a = sys.argv[1:]
out, title = a[0], a[1]
MM, NMIN, NMAX, MMS, NSMIN, NSMAX, NJ, NK, IVAC, NOWALL = (int(x) for x in a[2:12])
AWALL = float(a[12]); MODELK = int(a[13]); NSTA = int(a[14]); AL0 = float(a[15])
m_lo = int(a[16]) if len(a) > 16 else 1; m_hi = int(a[17]) if len(a) > 17 else MMS
MPINIT = int(a[18]) if len(a) > 18 else 0   # the stability table's column M stands for m = M + MPINIT
LPRESS = int(a[19]) if len(a) > 19 else 2   # 2: pressure gradient from the VMEC file; 0: reconstructed from the force balance
LCURRF = int(a[20]) if len(a) > 20 else 9   # 9: parallel current from the geometry; 1: from charge conservation, which loses the edge current of a pedestal (factor 8 in the growth rate of the cut cbm18a at n = 10)
L = []
L.append(' ' * 15 + '%-8s' % title); L.append('C'); L.append('C        MM  NMIN  NMAX   MMS NSMIN NSMAX NPROCS INSOL')
L.append('     ' + ''.join('%6d' % v for v in (MM, NMIN, NMAX, MMS, NSMIN, NSMAX, 1, 0)))
L.append('C'); L.append('C        NJ    NK  IVAC  LSSL  LSSD MMAXDF NMAXDF')
NMAXDF = max(38, 2 * max(abs(NSMIN), abs(NSMAX)) + 2)   # the convolution table holds sums of two toroidal numbers
L.append('     ' + ''.join('%6d' % v for v in (NJ, NK, IVAC, 4920, 1065)) + ''.join('%7d' % v for v in (230, NMAXDF)))
L += ['C', 'C     TABLE OF FOURIER COEFFIENTS FOR BOOZER COORDINATES', 'C     EQUILIBRIUM SETTINGS ARE COMPUTED FROM FIT/VMEC', 'C',
      'C M=  0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6  N']
for N in range(NMIN, NMAX + 1):
    row = [1 if m <= MM else 0 for m in range(37)]
    if N < 0: row[0] = 0
    L.append('     ' + ''.join('%2d' % v for v in row) + '%3d' % N)
L += ['C', '      LLAMPR      LVMTPR      LMETPR      LFOUPR', '           0           0           0           0',
      '      LLHSPR      LRHSPR      LEIGPR      LEFCPR', '           0           0           1           1',
      '      LXYZPR      LIOTPL      LDW2PL      LEFCPL', '           0           1           1           1',
      '      LCURRF      LMESHP      LMESHV      LITERS', '%12d' % LCURRF + '           1           2           1',
      '      LXYZPL      LEFPLS      LEQVPL      LPRESS', '           1           1           0' + '%12d' % LPRESS, 'C',
      'C    PVAC        PARFAC      QONAX        QN         DSVAC       QVAC    NOWALL']
L.append(''.join('%12.4E' % v for v in (1.6001, 0.0, 0.65, 0.0, 1.0, 1.6001)) + '     ' + '%2d' % NOWALL)
L += ['', 'C    AWALL       EWALL       DWALL       GWALL       DRWAL       DZWAL   NPWALL']
L.append(''.join('%12.4E' % v for v in (AWALL, 1.5, -1.5, 5.5, 0.0, 0.0)) + '     ' + '%2d' % 0)
L += ['C', 'C    RPLMIN       XPLO      DELTAJP      WCT        CURFAC']
L.append(''.join('%12.4E' % v for v in (1.0e-5, 1.0e-6, 4.0e-2, 0.66667, 1.0)))
L += ['C', 'C' + ' ' * 61 + 'MODELK = ' + '%6d' % MODELK, 'C',
      'C     NUMBER OF EQUILIBRIUM FIELD PERIODS PER STABILITY PERIOD: NSTA = ' + '%6d' % NSTA,
      'C', 'C     TABLE OF FOURIER COEFFIENTS FOR STABILITY DISPLACEMENTS', 'C',
      'C M=  0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5  N']
for N in range(NSMIN, NSMAX + 1):
    row = [1 if m_lo <= m <= m_hi else 0 for m in range(56)]
    L.append('     ' + ''.join('%2d' % v for v in row) + '%3d' % N)
L += ['C', 'C   NEV NITMAX         AL0      EPSPAM IGREEN MPINIT']
L.append('%7d%7d' % (1, 1500) + '%12.3E%12.4E' % (AL0, 1.0e-4) + '%7d%7d' % (0, MPINIT))
L.append('C')
open(out, 'w').write('\n'.join(L) + '\n'); print('wrote', out, len(L), 'lines')
