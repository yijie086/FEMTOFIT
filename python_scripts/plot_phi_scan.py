import numpy as np
import matplotlib.pyplot as plt
import gepard as g
from gepard.fits import th_KM15

# ==== 设置物理参数 ====
xB = 0.17
Q2 = 1.35
t = -0.21
E = 7.546
charge = -1
pol = -1

# ==== 扫描 φ（避免闭环）====
npoints = 90
phis = np.linspace(0, 2 * np.pi, npoints, endpoint=False)
phi_deg_shifted = (np.degrees(phis) + 180) % 360

xsecs = []
bsas = []

for phi in phis:
    pt = g.DataPoint(
        xB=xB,
        t=t,
        Q2=Q2,
        phi=phi,
        process='ep2epgamma',
        exptype='fixed target',
        in1energy=E,
        in1charge=charge,
        in1polarization=pol,
        observable='XS'
    )
    ptU = g.DataPoint(
        xB=xB,
        t=t,
        Q2=Q2,
        phi=phi,
        process='ep2epgamma',
        exptype='fixed target',
        in1energy=E,
        in1charge=charge,
        in1polarization=0,
        observable='XS'
    )
    pt.prepare()

    xsecs.append(th_KM15.XS(ptU))
    bsas.append(th_KM15.ALUI(pt))

# ==== Cross Section: 散点图 ====
plt.figure(figsize=(8, 6))
plt.scatter(phi_deg_shifted, xsecs, color='blue', label="KM15 XS", s=25)
plt.yscale("log")
plt.xlabel(r"$\phi$ [deg] (shifted)")
plt.ylabel(r"Cross Section [nb/GeV$^4$]")
plt.title(f"Cross Section vs shifted φ\nQ²={Q2}, xB={xB}, -t={-t}, E={E} GeV")
plt.grid(True)
plt.tight_layout()
plt.savefig("XS_phi_scatter_th_KM15.png", dpi=300)
plt.close()

# ==== BSA: 散点图 ====
plt.figure(figsize=(8, 6))
plt.scatter(phi_deg_shifted, bsas, color='red', label="KM15 BSA", s=25)
plt.xlabel(r"$\phi$ [deg] (shifted)")
plt.ylabel(r"Beam Spin Asymmetry (A$_{LU}$)")
plt.title(f"BSA vs shifted φ\nQ²={Q2}, xB={xB}, -t={-t}, E={E} GeV")
plt.grid(True)
plt.tight_layout()
plt.savefig("BSA_phi_scatter_th_KM15.png", dpi=300)
plt.close()

print("✅ 成功生成散点图：XS_phi_scatter_th_KM15.png 和 BSA_phi_scatter_th_KM15.png")
