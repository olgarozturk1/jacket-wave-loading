import numpy as np
import matplotlib.pyplot as plt

rho, CD, CM = 1025.0, 1.2, 0.8
D, a, d = 1.0, 6.0, 40.0
k, omega = 0.043, 0.628
length_factor = 1.23

# Drag: -35 m to +6 m
z = np.linspace(-35.0, 6.0, 4001)
weights = np.full(z.size, z[1] - z[0])
weights[[0, -1]] *= 0.5

# Inertia: -35 m to SWL (0 m)
z_inertia = np.linspace(-35.0, 0.0, 4001)
weights_inertia = np.full(
    z_inertia.size, z_inertia[1] - z_inertia[0]
)
weights_inertia[[0, -1]] *= 0.5

B = a * omega / np.sinh(k * d)
U = B * np.cosh(k * (z + d))
U_inertia = B * np.cosh(k * (z_inertia + d))

phase = np.linspace(0, 2 * np.pi, 8001)
time = phase / omega

# Independent-component horizontal drag: unchanged
Fd = np.zeros_like(phase)
for j, theta in enumerate(phase):
    u = U * np.sin(theta)
    Fd[j] = (
        length_factor * 0.5 * rho * CD * D
        * np.sum(weights * u * np.abs(u))
    )

# Horizontal inertia: upper limit changed to 0 m
Fi = (
    length_factor * rho * CM * np.pi * D**2 / 4
    * omega * np.sum(weights_inertia * U_inertia)
    * np.cos(phase)
)

# Add at the same phase
Ftotal = Fd + Fi

print("\nYELLOW BRACE: HORIZONTAL x FORCES FROM u")
print("Drag limits: -35 to +6 m")
print("Inertia limits: -35 to 0 m")
print("dl/dz = 1.23; x = 0")
print("\nForce       Maximum positive (kN)    Phase (degrees)")

for name, force in [
    ("Drag", Fd), ("Inertia", Fi), ("Total", Ftotal)
]:
    j = np.argmax(force)
    print(
        f"{name:<12}"
        f"{force[j] / 1000:>18.2f}"
        f"{np.degrees(phase[j]):>20.2f}"
    )

j = np.argmax(Ftotal)
print(
    f"\nAt total maximum: drag={Fd[j] / 1000:.2f} kN, "
    f"inertia={Fi[j] / 1000:.2f} kN, "
    f"time={time[j]:.3f} s"
)

plt.figure(figsize=(10, 5))

for name, force in [
    ("Drag", Fd), ("Inertia", Fi), ("Total", Ftotal)
]:
    plt.plot(np.degrees(phase), force / 1000, label=name)

plt.scatter(
    np.degrees(phase[j]),
    Ftotal[j] / 1000,
    color="red",
    zorder=5
)
plt.axhline(0, color="black", linewidth=0.7)
plt.xlabel("Wave phase at x=0 (degrees)")
plt.ylabel("Horizontal force (kN)")
plt.title("Yellow y-z brace: drag to crest; inertia to SWL")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()

