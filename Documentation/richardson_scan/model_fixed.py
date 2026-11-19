"""model_fixed.py <folder with onemode.py> -- the one-mode model of the theory document for the orders 1 to 3
at a fixed and at a free edge, with equal weights and with the weights of the sum rule."""
import sys, math
sys.path.insert(0, sys.argv[1])               # the folder "model" of the theory document, with onemode.py
import onemode as om
ns = (101, 201, 401, 801, 1601)
for edge, V, lab in (('fixed', om.V_edge, 'V_edge'), ('free', om.V_edge, 'V_edge'), ('fixed', om.V_flat, 'V_flat')):
    print(f"edge = {edge}, potential {lab}, true eigenvalue {om.reference(edge, V)[0]:.6f}")
    for p in (1, 2, 3):
        for kind in ('all 1', 'sum rule'):
            try:
                e = [float(x) for x in om.errors(ns, edge=edge, V=V, p=p, kind=kind)]
            except Exception as ex:
                print(f"  p={p} {kind}: {ex}"); continue
            orders = [math.log2(abs(e[i]/e[i+1])) if e[i+1] != 0 else float('nan') for i in range(len(e)-1)]
            print(f"  p={p} {kind:9s} errors " + " ".join(f"{x:+.2e}" for x in e) + "  orders " + " ".join(f"{o:.2f}" for o in orders) + f"   N*err(last) = {e[-1]*ns[-1]:+.3f}")
