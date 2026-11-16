import numpy as np
import matplotlib.pyplot as plt

# PARAMETERS
rho = 1025.0
CD = 1.2
CM = 0.8

D = 2.0             # Beam diameter (m)
L = 30.0            # Beam length along y (m)
z_beam = -35.0      # Elevation relative to SWL (m)
x_beam = 0.0        # Constant x coordinate (m)

a = 6.0             # Wave amplitude (m)
d = 40.0            # Water depth (m)
k = 0.043           # Wave number (1/m)
omega = 0.628       # Angular frequency (rad/s)

# ONE WAVE PERIOD
T = 2 * np.pi / omega
time = np.linspace(0, T, 20001)
theta = omega * time - k * x_beam

# VERTICAL VELOCITY AND ACCELERATION
W = a * omega * np.sinh(k * (z_beam + d)) / np.sinh(k * d)

w = W * np.cos(theta)
az = -omega * W * np.sin(theta)

# VERTICAL FORCES — LECTURE-NOTE FORMULATION
# Conditions are uniform along y, so integration gives a factor L.
Fd_z = 0.5 * rho * CD * D * L * w * np.abs(w)

Fi_z = rho * CM * (np.pi * D**2 / 4) * L * az

# ADD FORCES AT THE SAME TIME
Ftotal_z = Fd_z + Fi_z

# RESULTS
print("\nBEAM 1 ALONG y — VERTICAL FORCES")
print(f"Length: {L:.1f} m | Diameter: {D:.1f} m")
print(f"Elevation: z = {z_beam:.1f} m")
print("Lecture method: vertical drag uses w * abs(w).")
print("Positive force = upward; negative force = downward.")

for heading, upward in [
    ("MAXIMUM UPWARD FORCES", True),
    ("MAXIMUM DOWNWARD FORCES", False)
]:
    print(f"\n{heading}")
    print(f"{'Force':<12}{'Force (kN)':>14}"
          f"{'Time (s)':>12}{'Phase (deg)':>16}")

    for name, force in [
        ("Drag", Fd_z),
        ("Inertia", Fi_z),
        ("Total", Ftotal_z)
    ]:
        j = np.argmax(force) if upward else np.argmin(force)
        phase = np.degrees(theta[j]) % 360

        print(
            f"{name:<12}{force[j]/1000:>14.2f}"
            f"{time[j]:>12.3f}{phase:>16.2f}"
        )

# CONTRIBUTIONS AT MAXIMUM UPWARD TOTAL
idx_up = np.argmax(Ftotal_z)
idx_down = np.argmin(Ftotal_z)

print("\nAT MAXIMUM UPWARD TOTAL:")
print(f"Drag:    {Fd_z[idx_up]/1000:+.2f} kN")
print(f"Inertia: {Fi_z[idx_up]/1000:+.2f} kN")
print(f"Total:   {Ftotal_z[idx_up]/1000:+.2f} kN")

# PLOT
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
    color="red", zorder=5, label="Total extrema"
)

plt.axhline(0, color="black", linewidth=0.7)
plt.xlabel("Time (s)")
plt.ylabel("Vertical force (kN), positive upward")
plt.title("Beam 1 along y ")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()