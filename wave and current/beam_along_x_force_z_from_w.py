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

a = 6.0
d = 40.0
omega = 0.628       # Absolute angular frequency (rad/s)
U_current = 1.6     # Uniform following current along +x (m/s)

# Beam runs along x from 0 to 30 m.
# Current is parallel to the beam axis.
# Only transverse vertical Morison forces are calculated.
# Axial friction and end forces are not included.


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

# ---------- SPATIAL INTEGRATION ----------
x = np.linspace(0.0, L, 4001)

weights = np.full(x.size, x[1] - x[0])
weights[0] *= 0.5
weights[-1] *= 0.5

# ---------- ONE OBSERVED WAVE PERIOD ----------
T = 2 * np.pi / omega
time = np.linspace(0.0, T, 8001)
phase = omega * time

# Vertical wave velocity amplitude at beam elevation
W = (
    a * sigma
    * np.sinh(k * (z_beam + d))
    / np.sinh(k * d)
)

Fd_z = np.zeros_like(time)
Fi_z = np.zeros_like(time)

# ---------- FORCES ----------
for j, phi in enumerate(phase):
    theta = phi - k * x

    # Current has no direct vertical velocity component
    w = W * np.cos(theta)

    # Linearized material acceleration:
    # az = dw/dt + U_current * dw/dx
    #    = -(omega - k*U_current) * W * sin(theta)
    az = -sigma * W * np.sin(theta)

    Fd_z[j] = (
        0.5 * rho * CD * D
        * np.sum(weights * w * np.abs(w))
    )

    Fi_z[j] = (
        rho * CM * (np.pi * D**2 / 4)
        * np.sum(weights * az)
    )

# Add forces at the SAME time
Ftotal_z = Fd_z + Fi_z

# ---------- RESULTS ----------
print("\nBEAM 2 ALONG x — WAVE-CURRENT VERTICAL FORCES")
print(f"Integration: x = 0 to {L:.1f} m, at z = {z_beam:.1f} m")
print(f"Uniform current:      {U_current:.3f} m/s along +x")
print(f"Absolute frequency:   {omega:.6f} rad/s")
print(f"Intrinsic frequency:  {sigma:.6f} rad/s")
print(f"Wavenumber:           {k:.6f} 1/m")
print(f"Wavelength:           {2 * np.pi / k:.3f} m")
print(f"Observed period:      {T:.3f} s")
print("Inertia uses linearized material acceleration.")
print("Positive force = upward; negative force = downward.")
print("Reported phase is omega*t at x=0.")

for heading, upward in [
    ("MAXIMUM UPWARD FORCES", True),
    ("MAXIMUM DOWNWARD FORCES", False)
]:
    print(f"\n{heading}")
    print(
        f"{'Force':<12}"
        f"{'Force (kN)':>14}"
        f"{'Phase (deg)':>16}"
        f"{'Time (s)':>12}"
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
            f"{np.degrees(phase[idx]) % 360:>16.2f}"
            f"{time[idx]:>12.3f}"
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
    print(f"Phase:   {np.degrees(phase[idx]) % 360:.2f} degrees")
    print(f"Time:    {time[idx]:.3f} s")

# ---------- PLOT ----------
plt.figure(figsize=(10, 5))

plt.plot(
    np.degrees(phase), Fd_z / 1000,
    label="Vertical drag"
)
plt.plot(
    np.degrees(phase), Fi_z / 1000,
    label="Vertical inertia"
)
plt.plot(
    np.degrees(phase), Ftotal_z / 1000,
    label="Vertical total", linewidth=2
)

plt.scatter(
    np.degrees(phase[idx_up]),
    Ftotal_z[idx_up] / 1000,
    color="red", zorder=5,
    label="Maximum upward total"
)
plt.scatter(
    np.degrees(phase[idx_down]),
    Ftotal_z[idx_down] / 1000,
    color="darkorange", zorder=5,
    label="Maximum downward total"
)

plt.axhline(0, color="black", linewidth=0.7)
plt.xlabel("Reference wave phase at x=0 (degrees)")
plt.ylabel("Vertical force (kN), positive upward")
plt.title(
    "Beam 2 along x — linear wave-current interaction\n"
    "Uniform following current = 1.6 m/s"
)
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()