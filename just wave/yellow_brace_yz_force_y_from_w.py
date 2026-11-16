import numpy as np
import matplotlib.pyplot as plt

# ---------- PARAMETERS ----------
rho, CD, CM = 1025.0, 1.2, 0.8
D, a, d = 1.0, 6.0, 40.0
k, omega = 0.043, 0.628

alpha = np.arctan2(42.0, 30.0)  # Angle from horizontal
dl_dz = 1.23
x_brace = 0.0

# ---------- FIXED INTEGRATION LIMITS ----------
# Drag from w: -35 m to SWL (0 m)
z_drag = np.linspace(-35.0, 0.0, 20001)
weights_drag = np.full(
    z_drag.size, z_drag[1] - z_drag[0]
)
weights_drag[[0, -1]] *= 0.5

# Inertia from vertical acceleration: -35 m to crest (+6 m)
z_inertia = np.linspace(-35.0, 6.0, 20001)
weights_inertia = np.full(
    z_inertia.size, z_inertia[1] - z_inertia[0]
)
weights_inertia[[0, -1]] *= 0.5

# ---------- VELOCITY AMPLITUDES ----------
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
# Drag: integrated to SWL
AD = (
    0.5 * rho * CD * D * dl_dz
    * np.sum(weights_drag * Wn_drag * np.abs(Wn_drag))
)

# Inertia: integrated to crest
AI = (
    rho * CM * (np.pi * D**2 / 4) * dl_dz * omega
    * np.sum(weights_inertia * Wn_inertia)
)


def forces(x, t):
    """Horizontal y forces caused by w: drag, inertia, total, in N."""
    theta = omega * np.asarray(t) - k * np.asarray(x)

    # w varies as cos(theta); vertical acceleration as -sin(theta)
    Fd_n = AD * np.cos(theta) * np.abs(np.cos(theta))
    Fi_n = -AI * np.sin(theta)

    # Front brace falls in +y:
    # positive normal points towards +y and +z.
    Fd_y = Fd_n * np.sin(alpha)
    Fi_y = Fi_n * np.sin(alpha)

    return Fd_y, Fi_y, Fd_y + Fi_y


def main():
    T = 2 * np.pi / omega
    time = np.linspace(0, T, 100001, endpoint=False)

    Fd, Fi, Ft = forces(x_brace, time)
    phase = np.mod(
        np.degrees(omega * time - k * x_brace), 360
    )

    # ---------- RESULTS ----------
    print("\nYELLOW BRACE: HORIZONTAL y FORCE CAUSED BY w")
    print("Positive direction: +y.")
    print("Drag integration:    z = -35 to 0 m (SWL)")
    print("Inertia integration: z = -35 to +6 m (crest)")
    print(f"dl/dz = {dl_dz:.2f}")
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

    # ---------- PLOT ----------
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
    plt.ylabel("Horizontal force Fy (kN)")
    plt.title(
        "Yellow brace: horizontal y force caused by w\n"
        "Fixed limits: drag to SWL; inertia to crest"
    )
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
