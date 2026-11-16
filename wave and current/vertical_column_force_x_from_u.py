import numpy as np
import matplotlib.pyplot as plt

# ---------- PARAMETERS ----------
rho = 1025.0
g = 9.81
CD = 1.2
CM = 0.8
D = 4.0

a = 6.0                 # Local wave amplitude (m)
d = 40.0
omega = 0.628           # Absolute angular frequency (rad/s)
U_current = 1.6         # Uniform following current, +x (m/s)
x_column = 0.0

# Model:
# Linear wave-current interaction.
# Inertia uses linearized material acceleration about uniform current.
# Prescribed fixed limits: drag to crest, inertia to SWL.


# ---------- WAVE-CURRENT DISPERSION ----------
def solve_wavenumber(omega, current, depth, gravity):
    """Positive-intrinsic-frequency branch for following current."""
    if current <= 0:
        raise ValueError("This solver assumes a positive following current.")

    # sigma > 0 requires 0 < k < omega/current.
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
T = 2 * np.pi / omega


# ---------- INTEGRATION GRIDS ----------
z_drag = np.linspace(-d, a, 20001)
z_inertia = np.linspace(-d, 0.0, 20001)


def integrate_trapezoidal(values, z):
    return np.sum(
        0.5 * (values[:-1] + values[1:]) * np.diff(z)
    )


def wave_velocity_amplitude(z):
    return (
        a * sigma
        * np.cosh(k * (z + d))
        / np.sinh(k * d)
    )


U_drag = wave_velocity_amplitude(z_drag)
U_inertia = wave_velocity_amplitude(z_inertia)

# Linearized material acceleration amplitude = sigma * U(z)
inertia_amplitude = (
    rho * CM * (np.pi * D**2 / 4)
    * integrate_trapezoidal(sigma * U_inertia, z_inertia)
)


# ---------- TIME HISTORY ----------
time = np.linspace(0, T, 20001)
theta = omega * time - k * x_column
phase = np.degrees(theta) % 360

Fd = np.zeros_like(time)

for i, theta_i in enumerate(theta):
    u_wave = U_drag * np.sin(theta_i)
    u_total = U_current + u_wave

    Fd[i] = (
        0.5 * rho * CD * D
        * integrate_trapezoidal(
            u_total * np.abs(u_total), z_drag
        )
    )

Fi = inertia_amplitude * np.cos(theta)
Ftotal = Fd + Fi


# ---------- RESULTS ----------
print("\nVERTICAL COLUMN — LINEAR WAVE-CURRENT INTERACTION")
print(f"Following current:    {U_current:.3f} m/s")
print(f"Absolute frequency:   {omega:.6f} rad/s")
print(f"Intrinsic frequency:  {sigma:.6f} rad/s")
print(f"Wavenumber:           {k:.6f} 1/m")
print(f"Wavelength:           {2 * np.pi / k:.3f} m")
print(f"Observed period:      {T:.3f} s")
print(f"Drag limits:          {-d:.1f} to {a:.1f} m")
print(f"Inertia limits:       {-d:.1f} to 0.0 m")
print("Inertia uses linearized material acceleration.")
print("Signed forces below are in kN.")

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
    print(f"Phase:   {phase[idx]:.2f} degrees")


# ---------- PLOT ----------
plt.figure(figsize=(10, 5))

plt.plot(time, Fd / 1000, label="Wave-current drag")
plt.plot(time, Fi / 1000, label="Wave-current inertia")
plt.plot(time, Ftotal / 1000, label="Total", linewidth=2)

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
    "Vertical column — linear wave-current interaction\n"
    "Fixed limits: drag to crest; inertia to SWL"
)
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()