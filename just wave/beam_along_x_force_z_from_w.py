import numpy as np
import matplotlib.pyplot as plt

# PARAMETERS
rho = 1025.0
CD = 1.2
CM = 0.8

D = 2.0             # Beam diameter (m)
L = 30.0            # Beam length along x (m)
z_beam = -35.0      # Elevation relative to SWL (m)

a = 6.0
d = 40.0
k = 0.043
omega = 0.628

# COORDINATES:
# Beam runs from x=0 to x=30 m.
# Waves travel along +x.
# Positive vertical force points upward (+z).

# SPATIAL INTEGRATION POINTS
x = np.linspace(0.0, L, 4001)

# Trapezoidal integration weights
weights = np.full(x.size, x[1] - x[0])
weights[0] *= 0.5
weights[-1] *= 0.5

# ONE WAVE PERIOD
T = 2 * np.pi / omega
time = np.linspace(0.0, T, 8001)
phase = omega * time

# VERTICAL VELOCITY AMPLITUDE AT THE BEAM
W = a * omega * np.sinh(k * (z_beam + d)) / np.sinh(k * d)

Fd_z = np.zeros_like(time)
Fi_z = np.zeros_like(time)

# INTEGRATE ALONG THE BEAM AT EACH TIME
for j, phi in enumerate(phase):
    theta = phi - k * x

    w = W * np.cos(theta)
    az = -omega * W * np.sin(theta)

    # Lecture-note formula: signed vertical drag
    Fd_z[j] = (
        0.5 * rho * CD * D
        * np.sum(weights * w * np.abs(w))
    )

    Fi_z[j] = (
        rho * CM * (np.pi * D**2 / 4)
        * np.sum(weights * az)
    )

# TOTAL AT THE SAME TIME
Ftotal_z = Fd_z + Fi_z

# RESULTS
print("\nBEAM 2 ALONG x — VERTICAL FORCES")
print("Integration: x = 0 to 30 m, at z = -35 m")
print("Positive force = upward; negative force = downward.")
print("Reported phase is omega*t at x=0.")

for heading, upward in [
    ("MAXIMUM UPWARD FORCES", True),
    ("MAXIMUM DOWNWARD FORCES", False)
]:
    print(f"\n{heading}")
    print(f"{'Force':<12}{'Force (kN)':>14}"
          f"{'Phase (deg)':>16}{'Time (s)':>12}")

    for name, force in [
        ("Drag", Fd_z),
        ("Inertia", Fi_z),
        ("Total", Ftotal_z)
    ]:
        j = np.argmax(force) if upward else np.argmin(force)

        print(
            f"{name:<12}{force[j]/1000:>14.2f}"
            f"{np.degrees(phase[j]):>16.2f}"
            f"{time[j]:>12.3f}"
        )

idx = np.argmax(Ftotal_z)

print("\nAT MAXIMUM UPWARD TOTAL:")
print(f"Drag:    {Fd_z[idx]/1000:+.2f} kN")
print(f"Inertia: {Fi_z[idx]/1000:+.2f} kN")
print(f"Total:   {Ftotal_z[idx]/1000:+.2f} kN")

# PLOT
plt.figure(figsize=(10, 5))
plt.plot(np.degrees(phase), Fd_z/1000, label="Vertical drag")
plt.plot(np.degrees(phase), Fi_z/1000, label="Vertical inertia")
plt.plot(
    np.degrees(phase), Ftotal_z/1000,
    label="Vertical total", linewidth=2
)

plt.scatter(
    np.degrees(phase[idx]), Ftotal_z[idx]/1000,
    color="red", zorder=5, label="Maximum upward total"
)

plt.axhline(0, color="black", linewidth=0.7)
plt.xlabel("Reference wave phase at x=0 (degrees)")
plt.ylabel("Vertical force (kN), positive upward")
plt.title("Beam 2 along x ")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()
