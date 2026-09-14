import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# ЧТЕНИЕ RCS
# ============================================================

filename = "results/rcs.csv"

data = np.genfromtxt(
    filename,
    delimiter=",",
    names=True
)

theta = data["theta_deg"]
phi = data["phi_deg"]
sigma = data["sigma_m2"]
sigma_dbsm = data["sigma_dbsm"]


# ============================================================
# RCS [m^2]
# ============================================================

plt.figure()

plt.plot(phi, sigma)

plt.xlabel(r"$\phi$, deg")
plt.ylabel(r"$\sigma$, m$^2$")
plt.title(r"Radar Cross Section, $\theta=90^\circ$")

plt.grid(True)

plt.tight_layout()

plt.savefig(
    "results/rcs_m2.png",
    dpi=200
)


# ============================================================
# RCS [dBsm]
# ============================================================

plt.figure()

plt.plot(phi, sigma_dbsm)

plt.xlabel(r"$\phi$, deg")
plt.ylabel(r"$\sigma$, dBsm")
plt.title(r"Radar Cross Section, $\theta=90^\circ$")

plt.grid(True)

plt.tight_layout()

plt.savefig(
    "results/rcs_dbsm.png",
    dpi=200
)


# ============================================================
# ПОЛЯРНАЯ ДИАГРАММА RCS [dBsm]
# ============================================================

phi_rad = np.deg2rad(phi)

fig = plt.figure()

ax = fig.add_subplot(
    111,
    projection="polar"
)

ax.plot(
    phi_rad,
    sigma_dbsm
)

ax.set_theta_zero_location("E")
ax.set_theta_direction(1)

ax.set_title(
    r"RCS, $\theta=90^\circ$ [dBsm]"
)

ax.grid(True)

plt.tight_layout()

plt.savefig(
    "results/rcs_polar_dbsm.png",
    dpi=200
)


# ============================================================
# ПОКАЗАТЬ ГРАФИКИ
# ============================================================

plt.show()
