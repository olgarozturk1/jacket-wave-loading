import numpy as np
import matplotlib.pyplot as plt

# ---------- PARAMETERS ----------
rho, CD, CM = 1025.0, 1.2, 0.8
D, a, d = 1.0, 6.0, 40.0
k, omega = 0.043, 0.628

alpha = np.arctan2(42.0, 30.0)
dl_dz = 1.23
x_brace = 0.0

# ---------- FIXED INTEGRATION LIMITS ----------
# Drag from w: -35 m to SWL
z_drag = np.linspace(-35.0, 0.0, 20001)
weights_drag = np.full(
    z_drag.size, z_drag[1] - z_drag[0]
)
weights_drag[[0, -1]] *= 0.5

# Inertia from vertical acceleration: -35 m to crest
z_inertia = np.linspace(-35.0, 6.0, 20001)
weights_inertia = np.full(
    z_inertia.size, z_inertia[1] - z_inertia[0]
)
weights_inertia[[0, -1]] *= 0.5

# ---------- NORMAL VELOCITY AMPLITUDES ----------
W_drag = (
    a * omega * np.sinh(k * (z_drag + d))
    / np.sinh(k * d)
)
Wn_drag = W_drag * np.cos(alpha)

W_inertia = (
    a * omega * np.sinh(k * (z_inertia + d))
    / np.sinh(k * d)
)
Wn_inertia = W_inertia * np.cos(alpha)

# ---------- NORMAL FORCE AMPLITUDES ----------
AD = (
    0.5 * rho * CD * D * dl_dz
    * np.sum(weights_drag * Wn_drag * np.abs(Wn_drag))
)

AI = (
    rho * CM * (np.pi * D**2 / 4) * dl_dz * omega
    * np.sum(weights_inertia * Wn_inertia)
)


def forces(x, t):
    """Vertical z forces caused by w: drag, inertia, total, in N."""
    theta = omega * np.asarray(t) - k * np.asarray(x)

    Fd_n = AD * np.cos(theta) * np.abs(np.cos(theta))
    Fi_n = -AI * np.sin(theta)

    # Resolve normal forces into vertical components
    Fd_z = Fd_n * np.cos(alpha)
    Fi_z = Fi_n * np.cos(alpha)

    return Fd_z, Fi_z, Fd_z + Fi_z


def main():
    T = 2 * np.pi / omega
    time = np.linspace(0, T, 100001, endpoint=False)

    Fd, Fi, Ft = forces(x_brace, time)
    phase = np.mod(
        np.degrees(omega * time - k * x_brace), 360
    )

    print("\nYELLOW BRACE: VERTICAL z FORCE CAUSED BY w")
    print("Positive direction: upward.")
    print("Drag integration:    -35 to 0 m (SWL)")
    print("Inertia integration: -35 to +6 m (crest)")
    print("Phase = omega*t - k*x, measured locally at the brace.")

    print("\nMAXIMUM POSITIVE FORCES — each row has its own peak")
    print(
        f'{"Force":<12}'
        f'{"Maximum (kN)":>16}'
        f'{"Time (s)":>12}'
        f'{"Phase (deg)":>16}'
    )

    for name, values in [
        ("Drag", Fd),
        ("Inertia", Fi),
        ("Total", Ft)
    ]:
        idx = np.argmax(values)
        print(
            f"{name:<12}"
            f"{values[idx] / 1000:>16.2f}"
            f"{time[idx]:>12.3f}"
            f"{phase[idx]:>16.2f}"
        )

    idx_total = np.argmax(Ft)

    print("\nAT MAXIMUM TOTAL — values at the SAME instant")
    print(f"Drag    = {Fd[idx_total] / 1000:+.2f} kN")
    print(f"Inertia = {Fi[idx_total] / 1000:+.2f} kN")
    print(f"Total   = {Ft[idx_total] / 1000:+.2f} kN")
    print(f"Time    = {time[idx_total]:.3f} s")
    print(f"Phase   = {phase[idx_total]:.2f} degrees")

    plt.figure(figsize=(10, 5))

    for name, values in [
        ("Drag from w", Fd),
        ("Inertia from w", Fi),
        ("Total from w", Ft)
    ]:
        plt.plot(time, values / 1000, label=name)

    plt.scatter(
        time[idx_total],
        Ft[idx_total] / 1000,
        color="red",
        zorder=5,
        label="Maximum total"
    )

    plt.axhline(0, color="black", linewidth=0.7)
    plt.xlabel("Time (s)")
    plt.ylabel("Vertical force Fz (kN), positive upward")
    plt.title(
        "Yellow brace: vertical z force caused by w\n"
        "Fixed limits: drag to SWL; inertia to crest"
    )
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()

