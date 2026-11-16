import numpy as np
import matplotlib.pyplot as plt

# ---------- PARAMETERS ----------
rho, CD, CM = 1025.0, 1.2, 0.8
g = 9.81
D, a, d = 1.0, 6.0, 40.0

omega = 0.628       # Absolute angular frequency (rad/s)
U_current = 1.6     # Uniform following current along +x (m/s)

alpha = np.arctan2(42.0, 30.0)
dl_dz = 1.23
x_brace = 0.0

# Separate-component model:
# Calculate y-directed forces from vertical wave motion.
# Front brace falls in +y.
# Drag to SWL; inertia to crest.


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

# ---------- VERTICAL VELOCITY AMPLITUDES ----------
W_drag = (
    a * sigma * np.sinh(k * (z_drag + d))
    / np.sinh(k * d)
)
Wn_drag = W_drag * np.cos(alpha)

W_inertia = (
    a * sigma * np.sinh(k * (z_inertia + d))
    / np.sinh(k * d)
)
Wn_inertia = W_inertia * np.cos(alpha)

# ---------- NORMAL FORCE AMPLITUDES ----------
# Independent-component drag from w, integrated to SWL
AD = (
    0.5 * rho * CD * D * dl_dz
    * np.sum(weights_drag * Wn_drag * np.abs(Wn_drag))
)

# Linearized material acceleration amplitude:
# az = -sigma * W(z) * sin(theta)
# Project normally and integrate to crest.
AI = (
    rho * CM * (np.pi * D**2 / 4) * dl_dz * sigma
    * np.sum(weights_inertia * Wn_inertia)
)


def forces(x, t):
    """Horizontal y forces caused by w: drag, inertia, total, in N."""
    theta = omega * np.asarray(t) - k * np.asarray(x)

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
    print("Linear wave-current interaction; separate-component drag.")
    print("Positive direction: +y.")
    print(f"Following current:   {U_current:.3f} m/s along +x")
    print(f"Absolute frequency:  {omega:.6f} rad/s")
    print(f"Intrinsic frequency: {sigma:.6f} rad/s")
    print(f"Wavenumber:          {k:.6f} 1/m")
    print(f"Wavelength:          {2 * np.pi / k:.3f} m")
    print("Drag integration:    z = -35 to 0 m (SWL)")
    print("Inertia integration: z = -35 to +6 m (crest)")
    print(f"dl/dz = {dl_dz:.2f}")
    print("Phase = omega*t - k*x, measured locally at the brace.")

    for heading, positive in [
        ("MAXIMUM POSITIVE FORCES", True),
        ("MINIMUM SIGNED FORCES", False)
    ]:
        print(f"\n{heading} — each row has its own peak")
        print(
            f'{"Force":<12}'
            f'{"Force (kN)":>16}'
            f'{"Time (s)":>12}'
            f'{"Phase (deg)":>16}'
        )

        for name, values in [
            ("Drag", Fd),
            ("Inertia", Fi),
            ("Total", Ft)
        ]:
            idx = (
                np.argmax(values) if positive
                else np.argmin(values)
            )

            print(
                f"{name:<12}"
                f"{values[idx] / 1000:>16.2f}"
                f"{time[idx]:>12.3f}"
                f"{phase[idx]:>16.2f}"
            )

    idx_max = np.argmax(Ft)
    idx_min = np.argmin(Ft)

    for title, idx in [
        ("MAXIMUM TOTAL", idx_max),
        ("MINIMUM TOTAL", idx_min)
    ]:
        print(f"\nAT {title} — SAME INSTANT")
        print(f"Drag    = {Fd[idx] / 1000:+.2f} kN")
        print(f"Inertia = {Fi[idx] / 1000:+.2f} kN")
        print(f"Total   = {Ft[idx] / 1000:+.2f} kN")
        print(f"Time    = {time[idx]:.3f} s")
        print(f"Phase   = {phase[idx]:.2f} degrees")

    # ---------- PLOT ----------
    plt.figure(figsize=(10, 5))

    for name, values in [
        ("Drag from w", Fd),
        ("Inertia from w", Fi),
        ("Total from w", Ft)
    ]:
        plt.plot(time, values / 1000, label=name)

    plt.scatter(
        time[[idx_max, idx_min]],
        Ft[[idx_max, idx_min]] / 1000,
        color="red",
        zorder=5,
        label="Total extrema"
    )

    plt.axhline(0, color="black", linewidth=0.7)
    plt.xlabel("Time (s)")
    plt.ylabel("Horizontal force Fy (kN), positive in +y")
    plt.title(
        "Yellow brace — wave-current y force caused by w\n"
        "Fixed limits: drag to SWL; inertia to crest"
    )
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()