#!/usr/bin/env python3
"""grad gamma of the ring at the centroid of the boundary points, for the boundaries of the study
(analytic torus R_0 = 3, a = 1, as in exact_exterior_2d_options): circle, D-shape (kappa 1.7, delta 0.3), beans
R = R_0 + a (cos t + beta cos 2t), Z = 1.5 a sin t with beta = 0.6 and 0.9.
gamma = 1 + |x - x_s|^2 / (2 R R_s),  grad gamma = (x - x_s)/(R R_s) - (gamma - 1) e_R / R.
Option A's row i is singular where the outward normal n is perpendicular to grad gamma.
Writes 4_grad_gamma.png.   Usage: python3 plot_grad_gamma.py <folder of this script>"""
import sys, os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

R0, a = 3.0, 1.0
def shape(name, t):
    if name == 'circle':
        return R0 + a*np.cos(t), a*np.sin(t)
    if name == 'D':
        return R0 + a*np.cos(t + 0.3*np.sin(t)), 1.7*a*np.sin(t)
    b = float(name.split('=')[1])
    return R0 + a*(np.cos(t) + b*np.cos(2*t)), 1.5*a*np.sin(t)

def gamma(R, Z, Rs, Zs):
    return 1 + ((R-Rs)**2 + (Z-Zs)**2)/(2*R*Rs)

def grad_gamma(R, Z, Rs, Zs):
    g = gamma(R, Z, Rs, Zs)
    return (R-Rs)/(R*Rs) - (g-1)/R, (Z-Zs)/(R*Rs)

cases = [('circle', 'circle'), ('D', 'D-shape (κ 1.7, δ 0.3)'), ('bean β=0.6', 'bean, β = 0.6'),
         ('bean β=0.9', 'bean, β = 0.9')]
cols = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100']
fig = plt.figure(figsize=(12, 7.2))
gs = fig.add_gridspec(2, 4, height_ratios=[1.35, 1])
t = np.linspace(0, 2*np.pi, 4001)
tp = np.linspace(0, 2*np.pi, 1201)[:-1]                 # points as PB3D would hold them
axl = fig.add_subplot(gs[1, :])
for i, (name, title) in enumerate(cases):
    ax = fig.add_subplot(gs[0, i])
    Rp, Zp = shape(name, tp); Rs, Zs = Rp.mean(), Zp.mean()   # ring at the centroid of the points
    Rb, Zb = shape(name, t)
    # contours of gamma: grad gamma is perpendicular to them
    Rg, Zg = np.meshgrid(np.linspace(1.6, 5.1, 400), np.linspace(-2.0, 2.0, 400))
    G = gamma(Rg, Zg, Rs, Zs)
    ax.contour(Rg, Zg, np.sqrt(G-1), levels=np.linspace(0.05, 0.9, 14), colors='#c9c9c9', linewidths=0.6)
    # direction of grad gamma on the boundary
    k = np.arange(0, len(t)-1, 160)
    gR, gZ = grad_gamma(Rb[k], Zb[k], Rs, Zs); gn = np.hypot(gR, gZ)
    ax.quiver(Rb[k], Zb[k], gR/gn, gZ/gn, color='#6a6a6a', scale=14, width=0.006, zorder=3)
    ax.plot(Rb, Zb, color=cols[i], lw=2.2, zorder=2)
    ax.plot(Rs, Zs, 'x', color='#222222', ms=8, mew=2, zorder=4)
    # outward normal and its angle with grad gamma
    dR, dZ = np.gradient(Rb, t), np.gradient(Zb, t)
    nR, nZ = dZ, -dR
    gR, gZ = grad_gamma(Rb, Zb, Rs, Zs)
    c = (nR*gR + nZ*gZ)/np.hypot(nR, nZ)/np.hypot(gR, gZ)
    z = np.where(np.diff(np.sign(c)) != 0)[0]
    ax.plot(Rb[z], Zb[z], 'o', ms=9, mfc='none', mec='#d03b3b', mew=2.2, zorder=5)
    ax.set_title(title, fontsize=10)
    ax.set_aspect('equal'); ax.set_xlim(1.6, 5.1); ax.set_ylim(-2.0, 2.0)
    ax.set_xlabel('R'); ax.tick_params(labelsize=8)
    if i == 0: ax.set_ylabel('Z')
    axl.plot(t/np.pi, c, color=cols[i], lw=2, label=title + '   (smallest: %+.2f)' % c.min())
axl.axhline(0, color='#d03b3b', lw=1)
axl.text(0.02, 0.04, 'A singular where the curve crosses 0', color='#d03b3b', fontsize=9)
axl.set_xlim(0, 2); axl.set_ylim(-1.05, 1.05)
axl.set_xlabel('t / π  (t = 0 outboard midplane, t = 1 inboard midplane)')
axl.set_ylabel('cos of the angle between n and ∇γ')
axl.grid(True, color='#e6e6e6', lw=0.6)
axl.legend(loc='lower left', fontsize=9, frameon=False)
fig.suptitle('∇γ of the ring at the centroid (×): arrows on the boundary, grey curves γ = const;'
             ' red circles where n is perpendicular to ∇γ', fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.97))
out = os.path.join(sys.argv[1], '4_grad_gamma.png')
fig.savefig(out, dpi=150); print(out)
