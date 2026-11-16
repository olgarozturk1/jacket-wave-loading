"""Beam along y: horizontal x forces from u.
Separate-component drag follows Morisons lecture, printed page 83.
"""
import numpy as np
import matplotlib.pyplot as plt

# PARAMETERS
rho = 1025.0
CD = 1.2
CM = 0.8

D = 2.0           # Horizontal member diameter (m)
L = 30.0          # Member length along y (m)
z_member = -35.0  # Member elevation relative to SWL (m)
x_member = 0.0    # Constant x position (m)

a = 6.0
d = 40.0         # Water depth (m)
k = 0.043
omega = 0.628

# TIME: ONE WAVE PERIOD
T = 2 * np.pi / omega
time = np.linspace(0, T, 20001)
theta = omega * time - k * x_member

# HORIZONTAL WAVE VELOCITY AT THE MEMBER
B = a * omega / np.sinh(k * d)
R = k * (z_member + d)

u = B * np.cosh(R) * np.sin(theta)

# HORIZONTAL ACCELERATION
ax = omega * B * np.cosh(R) * np.cos(theta)

# HORIZONTAL FORCES: SEPARATE u|u| DRAG MODEL
# Integration along y reduces to multiplication by L.
Fd = 0.5 * rho * CD * D * L * u * np.abs(u)

Fi = rho * CM * (np.pi * D**2 / 4) * L * ax

# TOTAL AT EACH INSTANT
Ftotal = Fd + Fi

# MAXIMUM POSITIVE TOTAL
idx = np.argmax(Ftotal)
t_max = time[idx]
phase_max = np.degrees(theta[idx]) % 360

# RESULTS
print("\nHORIZONTAL MEMBER ALONG y — HORIZONTAL FORCES")
print(f"Length: {L:.1f} m | Diameter: {D:.1f} m")
print(f"Elevation: z = {z_member:.1f} m")
print("All forces below are in kN.\n")

print(f"Maximum drag alone:     {np.max(Fd) / 1000:10.2f}")
print(f"Maximum inertia alone:  {np.max(Fi) / 1000:10.2f}")

print("\nAT THE MAXIMUM TOTAL FORCE:")
print(f"Drag contribution:      {Fd[idx] / 1000:10.2f}")
print(f"Inertia contribution:   {Fi[idx] / 1000:10.2f}")
print(f"TOTAL FORCE:            {Ftotal[idx] / 1000:10.2f}")
print(f"\nTime:                   {t_max:.3f} s")
print(f"Local phase:            {phase_max:.2f} degrees")

# PLOT
plt.figure(figsize=(10, 5))

plt.plot(time, Fd / 1000, label="Horizontal drag")
plt.plot(time, Fi / 1000, label="Horizontal inertia")
plt.plot(
    time, Ftotal / 1000,
    label="Horizontal total", linewidth=2
)

plt.scatter(
    t_max, Ftotal[idx] / 1000,
    color="red", zorder=5,
    label=f"Maximum total = {Ftotal[idx] / 1000:.2f} kN"
)

plt.axhline(0, color="black", linewidth=0.7)
plt.xlabel("Time (s)")
plt.ylabel("Horizontal force (kN)")
plt.title("Lower horizontal member along y: drag + inertia")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()