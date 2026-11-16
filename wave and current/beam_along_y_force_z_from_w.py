import numpy as np
import matplotlib.pyplot as plt

# ---------- PARAMETERS ----------
rho = 1025.0
g = 9.81
CD = 1.2
CM = 0.8

D = 2.0
L = 30.0
z_beam = -35.0
x_beam = 0.0

a = 6.0
d = 40.0
omega = 0.628       # Absolute angular frequency (rad/s)
U_current = 1.6     # Uniform following current along +x (m/s)

# Beam lies along y at constant x and z.
# Vertical drag retains the separate-component formula w*abs(w).
# Inertia uses linearized material acceleration about uniform current.


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
theta = omega * time - k * x_beam
phase = np.degrees(theta) % 360

# ---------- VERTICAL WAVE KINEMATICS ----------
W = (
    a * sigma
    * np.sinh(k * (z_beam + d))
    / np.sinh(k * d)
)

w = W * np.cos(theta)

# Linearized material acceleration:
# az = dw/dt + U_current*dw/dx = -sigma*W*sin(theta)
az = -sigma * W * np.sin(theta)

# ---------- VERTICAL FORCES ----------
# Conditions are uniform along y, so integration gives a factor L.
# Separate-component drag: retain the original w*abs(w) model.
Fd_z = 0.5 * rho * CD * D * L * w * np.abs(w)

Fi_z = rho * CM * (np.pi * D**2 / 4) * L * az

# Add forces at the SAME time
Ftotal_z = Fd_z + Fi_z

# ---------- RESULTS ----------
print("\nBEAM 1 ALONG y — WAVE-CURRENT VERTICAL FORCES")
print(f"Length: {L:.1f} m | Diameter: {D:.1f} m")
print(f"Elevation: z = {z_beam:.1f} m | x = {x_beam:.1f} m")
print(f"Uniform current:     {U_current:.3f} m/s along +x")
print(f"Absolute frequency:  {omega:.6f} rad/s")
print(f"Intrinsic frequency: {sigma:.6f} rad/s")
print(f"Wavenumber:          {k:.6f} 1/m")
print(f"Wavelength:          {2 * np.pi / k:.3f} m")
print("Vertical drag uses the separate-component formula w*abs(w).")
print("Inertia uses linearized material acceleration.")
print("Positive force = upward; negative force = downward.")

for heading, upward in [
    ("MAXIMUM UPWARD FORCES", True),
    ("MAXIMUM DOWNWARD FORCES", False)
]:
    print(f"\n{heading}")
    print(
        f"{'Force':<12}"
        f"{'Force (kN)':>14}"
        f"{'Time (s)':>12}"
        f"{'Phase (deg)':>16}"
    )

    for name, force in [
        ("Drag", Fd_z),
        ("Inertia", Fi_z),
        ("Total", Ftotal_z)
    ]:
        idx = np.argmax(force) if upward else np.argmin(force)

        print(
            f"{name:<12}"
            f"{force[idx] / 1000:>14.2f}"
            f"{time[idx]:>12.3f}"
            f"{phase[idx]:>16.2f}"
        )

idx_up = np.argmax(Ftotal_z)
idx_down = np.argmin(Ftotal_z)

for title, idx in [
    ("MAXIMUM UPWARD TOTAL", idx_up),
    ("MAXIMUM DOWNWARD TOTAL", idx_down)
]:
    print(f"\nAT {title} — SAME INSTANT")
    print(f"Drag:    {Fd_z[idx] / 1000:+.2f} kN")
    print(f"Inertia: {Fi_z[idx] / 1000:+.2f} kN")
    print(f"Total:   {Ftotal_z[idx] / 1000:+.2f} kN")
    print(f"Time:    {time[idx]:.3f} s")
    print(f"Phase:   {phase[idx]:.2f} degrees")

# ---------- PLOT ----------
plt.figure(figsize=(10, 5))

plt.plot(time, Fd_z / 1000, label="Vertical drag")
plt.plot(time, Fi_z / 1000, label="Vertical inertia")
plt.plot(
    time, Ftotal_z / 1000,
    label="Vertical total", linewidth=2
)

plt.scatter(
    time[[idx_up, idx_down]],
    Ftotal_z[[idx_up, idx_down]] / 1000,
    color="red", zorder=5,
    label="Total extrema"
)

plt.axhline(0, color="black", linewidth=0.7)
plt.xlabel("Time (s)")
plt.ylabel("Vertical force (kN), positive upward")
plt.title(
    "Beam 1 along y — linear wave-current interaction\n"
    "Vertical drag model: w|w|"
)
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()
