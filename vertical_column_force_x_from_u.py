import numpy as np
import matplotlib.pyplot as plt

# PARAMETERS
rho = 1025.0
CD = 1.2
CM = 0.8
D = 4.0
a = 6.0
d = 40.0
k = 0.043
omega = 0.628
x_column = 0.0

# SEPARATE FIXED INTEGRATION LIMITS
z_drag = np.linspace(-d, a, 20001)       # Seabed to crest
z_inertia = np.linspace(-d, 0.0, 20001)  # Seabed to SWL


def integrate_trapezoidal(values, z):
    return np.sum(
        0.5 * (values[:-1] + values[1:]) * np.diff(z)
    )


def velocity_amplitude(z):
    return a * omega * np.cosh(k * (z + d)) / np.sinh(k * d)


# HORIZONTAL VELOCITY AMPLITUDES
U_drag = velocity_amplitude(z_drag)
U_inertia = velocity_amplitude(z_inertia)

# INTEGRATED FORCE AMPLITUDES
drag_amplitude = (
    0.5 * rho * CD * D
    * integrate_trapezoidal(U_drag**2, z_drag)
)

inertia_amplitude = (
    rho * CM * (np.pi * D**2 / 4)
    * omega
    * integrate_trapezoidal(U_inertia, z_inertia)
)

# FORCES OVER ONE WAVE PERIOD
T = 2 * np.pi / omega
time = np.linspace(0, T, 20001)
theta = omega * time - k * x_column

# u = U sin(theta), du/dt = omega U cos(theta)
Fd = drag_amplitude * np.sin(theta) * np.abs(np.sin(theta))
Fi = inertia_amplitude * np.cos(theta)

# ADD FORCES AT THE SAME TIME
Ftotal = Fd + Fi

# MAXIMUM POSITIVE TOTAL FORCE
idx = np.argmax(Ftotal)
t_max = time[idx]
phase_max = np.degrees(theta[idx]) % 360

# RESULTS
print("\nVERTICAL COLUMN — HORIZONTAL FORCES")
print(f"Drag integration:    {-d:.1f} m to {a:.1f} m (crest)")
print(f"Inertia integration: {-d:.1f} m to 0.0 m (SWL)")
print("Assumption: separate fixed integration limits.")
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
plt.plot(time, Fd / 1000, label="Horizontal drag (to crest)")
plt.plot(time, Fi / 1000, label="Horizontal inertia (to SWL)")
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
plt.title("Vertical column — separate fixed integration limits")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()
