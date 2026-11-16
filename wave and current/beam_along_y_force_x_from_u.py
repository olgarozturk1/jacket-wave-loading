"""
Beam along y: horizontal x forces.
Linear wave-current interaction with uniform following current.
Separate-component drag: u_total * abs(u_total).
Inertia: linearized material acceleration about uniform current.
"""
import numpy as np
import matplotlib.pyplot as plt

# ---------- PARAMETERS ----------
rho = 1025.0
g = 9.81
CD = 1.2
CM = 0.8

D = 2.0
L = 30.0
z_member = -35.0
x_member = 0.0

a = 6.0
d = 40.0
omega = 0.628       # Absolute angular frequency (rad/s)
U_current = 1.6     # Uniform following current along +x (m/s)


# ---------- WAVE-CURRENT DISPERSION ----------
def solve_wavenumber(omega, current, depth, gravity):
    """Positive-intrinsic-frequency branch for following current."""
    if current <= 0:
        raise ValueError(
            "This solver assumes a positive following current."
        )

    lo = 0.0
    hi = omega / current

    for _ in range(100):
        mid = 0.5 * (lo + hi)

        residual = (
            (omega - current * mid)**2
            - gravity * mid * np.tanh(mid * depth)
        )

        if residual > 0:
            lo = mid
        else:
            hi = mid

    return 0.5 * (lo + hi)


k = solve_wavenumber(omega, U_current, d, g)
sigma = omega - k * U_current

# ---------- ONE OBSERVED WAVE PERIOD ----------
T = 2 * np.pi / omega
time = np.linspace(0, T, 20001)

theta = omega * time - k * x_member
phase = np.degrees(theta) % 360

# ---------- HORIZONTAL VELOCITY ----------
Au = (
    a * sigma
    * np.cosh(k * (z_member + d))
    / np.sinh(k * d)
)

u_wave = Au * np.sin(theta)
u_total = U_current + u_wave

# ---------- HORIZONTAL ACCELERATION ----------
# Linearized material acceleration:
# ax = du_wave/dt + U_current * du_wave/dx
#    = (omega - k*U_current) * Au * cos(theta)
ax = sigma * Au * np.cos(theta)

# ---------- HORIZONTAL FORCES ----------
# Uniform conditions along y: integration gives a factor L.
Fd = (
    0.5 * rho * CD * D * L
    * u_total * np.abs(u_total)
)

Fi = rho * CM * (np.pi * D**2 / 4) * L * ax

# Add forces at the SAME instant
Ftotal = Fd + Fi

idx_max = np.argmax(Ftotal)
idx_min = np.argmin(Ftotal)

# ---------- RESULTS ----------
print("\nBEAM ALONG y — HORIZONTAL x FORCES")
print("Linear wave-current interaction; separate-component drag.")
print(f"Length: {L:.1f} m | Diameter: {D:.1f} m")
print(f"Elevation: z = {z_member:.1f} m")
print(f"Current:             {U_current:.3f} m/s along +x")
print(f"Absolute frequency:  {omega:.6f} rad/s")
print(f"Intrinsic frequency: {sigma:.6f} rad/s")
print(f"Wavenumber:          {k:.6f} 1/m")
print(f"Wavelength:          {2 * np.pi / k:.3f} m")
print(f"Observed period:     {T:.3f} s")
print("All forces below are signed, in kN.")

print(
    f"\n{'Force':<12}"
    f"{'Max (kN)':>14}"
    f"{'Phase (deg)':>15}"
    f"{'Min (kN)':>14}"
    f"{'Phase (deg)':>15}"
)

for name, values in [
    ("Drag", Fd),
    ("Inertia", Fi),
    ("Total", Ftotal)
]:
    i_max = np.argmax(values)
    i_min = np.argmin(values)

    print(
        f"{name:<12}"
        f"{values[i_max] / 1000:>14.2f}"
        f"{phase[i_max]:>15.2f}"
        f"{values[i_min] / 1000:>14.2f}"
        f"{phase[i_min]:>15.2f}"
    )

for title, idx in [
    ("MAXIMUM TOTAL FORCE", idx_max),
    ("MINIMUM TOTAL FORCE", idx_min)
]:
    print(f"\nAT {title} — SAME INSTANT")
    print(f"Drag contribution:    {Fd[idx] / 1000:+.2f} kN")
    print(f"Inertia contribution: {Fi[idx] / 1000:+.2f} kN")
    print(f"TOTAL FORCE:          {Ftotal[idx] / 1000:+.2f} kN")
    print(f"Time:                 {time[idx]:.3f} s")
    print(f"Local phase:          {phase[idx]:.2f} degrees")

# ---------- PLOT ----------
plt.figure(figsize=(10, 5))

plt.plot(time, Fd / 1000, label="Wave-current horizontal drag")
plt.plot(time, Fi / 1000, label="Horizontal inertia")
plt.plot(
    time, Ftotal / 1000,
    label="Horizontal total", linewidth=2
)

plt.scatter(
    time[idx_max], Ftotal[idx_max] / 1000,
    color="red", zorder=5,
    label=f"Maximum total = {Ftotal[idx_max] / 1000:.2f} kN"
)
plt.scatter(
    time[idx_min], Ftotal[idx_min] / 1000,
    color="darkorange", zorder=5,
    label=f"Minimum total = {Ftotal[idx_min] / 1000:.2f} kN"
)

plt.axhline(0, color="black", linewidth=0.7)
plt.xlabel("Time (s)")
plt.ylabel("Horizontal force (kN), positive in +x")
plt.title(
    "Beam along y — horizontal x forces\n"
    "Linear wave-current interaction, current = 1.6 m/s"
)
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()