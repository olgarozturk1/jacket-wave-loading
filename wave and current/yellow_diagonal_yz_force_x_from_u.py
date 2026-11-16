import numpy as np
import matplotlib.pyplot as plt

# ---------- PARAMETERS ----------
rho, CD, CM = 1025.0, 1.2, 0.8
g = 9.81
D, a, d = 1.0, 6.0, 40.0

omega = 0.628       # Absolute angular frequency (rad/s)
U_current = 1.6     # Uniform following current along +x (m/s)
length_factor = 1.23

# Yellow brace lies in the y-z plane at x=0.
# Horizontal flow is entirely normal to the member.
# Separate-component horizontal drag is retained.


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

# ---------- FIXED INTEGRATION LIMITS ----------
# Drag: -35 m to crest (+6 m)
z_drag = np.linspace(-35.0, 6.0, 4001)
weights_drag = np.full(
    z_drag.size, z_drag[1] - z_drag[0]
)
weights_drag[[0, -1]] *= 0.5

# Inertia: -35 m to SWL (0 m)
z_inertia = np.linspace(-35.0, 0.0, 4001)
weights_inertia = np.full(
    z_inertia.size, z_inertia[1] - z_inertia[0]
)
weights_inertia[[0, -1]] *= 0.5

# ---------- WAVE VELOCITY AMPLITUDES ----------
B = a * sigma / np.sinh(k * d)

U_drag = B * np.cosh(k * (z_drag + d))
U_inertia = B * np.cosh(k * (z_inertia + d))

# At x=0, local phase is omega*t everywhere along the brace
phase = np.linspace(0, 2 * np.pi, 8001)
time = phase / omega

# ---------- HORIZONTAL DRAG ----------
Fd = np.zeros_like(phase)

for j, theta in enumerate(phase):
    u_wave = U_drag * np.sin(theta)
    u_total = U_current + u_wave

    Fd[j] = (
        length_factor * 0.5 * rho * CD * D
        * np.sum(
            weights_drag * u_total * np.abs(u_total)
        )
    )

# ---------- HORIZONTAL INERTIA ----------
# Linearized material acceleration:
# ax = (partial_t + U_current*partial_x) u_wave
#    = sigma * U(z) * cos(theta)
Fi = (
    length_factor * rho * CM * (np.pi * D**2 / 4)
    * sigma * np.sum(weights_inertia * U_inertia)
    * np.cos(phase)
)

# Add forces at the SAME phase
Ftotal = Fd + Fi

# ---------- RESULTS ----------
print("\nYELLOW y-z BRACE: HORIZONTAL x FORCES")
print("Linear wave-current interaction; separate-component drag.")
print(f"Uniform current:     {U_current:.3f} m/s along +x")
print(f"Absolute frequency:  {omega:.6f} rad/s")
print(f"Intrinsic frequency: {sigma:.6f} rad/s")
print(f"Wavenumber:          {k:.6f} 1/m")
print(f"Wavelength:          {2 * np.pi / k:.3f} m")
print("Drag limits: -35 to +6 m")
print("Inertia limits: -35 to 0 m")
print("dl/dz = 1.23; x = 0")
print("All forces below are signed, in kN.")

print(
    f"\n{'Force':<12}"
    f"{'Max (kN)':>14}"
    f"{'Phase (deg)':>15}"
    f"{'Min (kN)':>14}"
    f"{'Phase (deg)':>15}"
)

for name, force in [
    ("Drag", Fd),
    ("Inertia", Fi),
    ("Total", Ftotal)
]:
    i_max = np.argmax(force)
    i_min = np.argmin(force)

    print(
        f"{name:<12}"
        f"{force[i_max] / 1000:>14.2f}"
        f"{np.degrees(phase[i_max]) % 360:>15.2f}"
        f"{force[i_min] / 1000:>14.2f}"
        f"{np.degrees(phase[i_min]) % 360:>15.2f}"
    )

idx_max = np.argmax(Ftotal)
idx_min = np.argmin(Ftotal)

for title, idx in [
    ("MAXIMUM TOTAL", idx_max),
    ("MINIMUM TOTAL", idx_min)
]:
    print(f"\nAT {title} — SAME INSTANT")
    print(f"Drag:    {Fd[idx] / 1000:+.2f} kN")
    print(f"Inertia: {Fi[idx] / 1000:+.2f} kN")
    print(f"Total:   {Ftotal[idx] / 1000:+.2f} kN")
    print(f"Time:    {time[idx]:.3f} s")
    print(f"Phase:   {np.degrees(phase[idx]) % 360:.2f} degrees")

# ---------- PLOT ----------
plt.figure(figsize=(10, 5))

for name, force in [
    ("Drag", Fd),
    ("Inertia", Fi),
    ("Total", Ftotal)
]:
    plt.plot(
        np.degrees(phase),
        force / 1000,
        label=name,
        linewidth=2 if name == "Total" else 1.5
    )

plt.scatter(
    np.degrees(phase[[idx_max, idx_min]]),
    Ftotal[[idx_max, idx_min]] / 1000,
    color="red",
    zorder=5,
    label="Total extrema"
)

plt.axhline(0, color="black", linewidth=0.7)
plt.xlabel("Wave phase at x=0 (degrees)")
plt.ylabel("Horizontal force (kN), positive in +x")
plt.title(
    "Yellow y-z brace — wave-current horizontal forces\n"
    "Fixed limits: drag to crest; inertia to SWL"
)
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()