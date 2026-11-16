import numpy as np
import matplotlib.pyplot as plt

# ---------- PARAMETERS ----------
rho = 1025.0
g = 9.81
Cd = 1.2
Cm = 0.8
D = 1.0

a = 6.0
omega = 0.628       # Absolute angular frequency (rad/s)
d = 40.0
U_current = 1.6     # Following current along +x (m/s)

z_min = -35.0
z_swl = 0.0
z_crest = 6.0

rise = 42.0
span = 30.0

alpha = np.arctan2(rise, span)
s = np.sin(alpha)
c = np.cos(alpha)
m = span / rise

# Origin: front-left column at SWL.
# Brace rises from (x, z) = (0, -35) to (30, +7) m.
# Therefore x = m*(z + 35); its SWL crossing is at x = 25 m.
# Normal unit vector: n = (s, 0, -c).
#
# Preserved assumptions:
# Drag: complete normal velocity over -35 to +6 m.
# Inertia: ax to SWL; az to crest.


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
T = 2 * np.pi / omega


# ---------- DEPTH AND TIME GRIDS ----------
z_lower = np.linspace(z_min, z_swl, 4001)
z_upper = np.linspace(z_swl, z_crest, 2001)

# Global x-coordinates along the brace; depth grids are unchanged.
x_lower = m * (z_lower + 35.0)
x_upper = m * (z_upper + 35.0)
t_list = np.linspace(0, T, 8000, endpoint=False)


def integrate(values, z):
    return np.sum(
        0.5 * (values[:-1] + values[1:]) * np.diff(z)
    )


B = a * sigma / np.sinh(k * d)

U_lower = B * np.cosh(k * (z_lower + d))
W_lower = B * np.sinh(k * (z_lower + d))

U_upper = B * np.cosh(k * (z_upper + d))
W_upper = B * np.sinh(k * (z_upper + d))

Fd_x = np.zeros_like(t_list)
Fi_x = np.zeros_like(t_list)
Fd_z = np.zeros_like(t_list)
Fi_z = np.zeros_like(t_list)


# ---------- CALCULATE FORCES ----------
for i, t in enumerate(t_list):

    # ----- BELOW SWL -----
    theta_lower = omega * t - k * x_lower

    u_lower = U_current + U_lower * np.sin(theta_lower)
    w_lower = W_lower * np.cos(theta_lower)

    # Linearized material accelerations:
    # (partial_t + U_current * partial_x) of wave velocity
    ax_lower = sigma * U_lower * np.cos(theta_lower)
    az_lower = -sigma * W_lower * np.sin(theta_lower)

    un_lower = u_lower * s - w_lower * c
    an_lower = ax_lower * s - az_lower * c

    # ----- ABOVE SWL -----
    theta_upper = omega * t - k * x_upper

    u_upper = U_current + U_upper * np.sin(theta_upper)
    w_upper = W_upper * np.cos(theta_upper)

    # Both wave velocity components retained up to crest
    un_upper = u_upper * s - w_upper * c

    # Preserve the prescribed inertia cutoff:
    # ax excluded above SWL, az retained to crest
    az_upper = -sigma * W_upper * np.sin(theta_upper)
    an_upper = -az_upper * c

    # ----- NORMAL FORCE INTEGRALS -----
    I_drag = (
        integrate(un_lower * np.abs(un_lower), z_lower)
        + integrate(un_upper * np.abs(un_upper), z_upper)
    )

    I_inertia = (
        integrate(an_lower, z_lower)
        + integrate(an_upper, z_upper)
    )

    # dl = dz / sin(alpha)
    Fd_normal = 0.5 * rho * Cd * D * I_drag / s
    Fi_normal = rho * Cm * (np.pi * D**2 / 4) * I_inertia / s

    # Resolve into global components
    Fd_x[i] = Fd_normal * s
    Fi_x[i] = Fi_normal * s

    Fd_z[i] = -Fd_normal * c
    Fi_z[i] = -Fi_normal * c


# ---------- TOTAL AT THE SAME INSTANT ----------
Ft_x = Fd_x + Fi_x
Ft_z = Fd_z + Fi_z


# ---------- RESULTS ----------
def print_results(title, force_arrays):
    print("\n" + title)
    print(
        f"{'Force':<12}"
        f"{'Max (kN)':>14}"
        f"{'Phase (deg)':>15}"
        f"{'Min (kN)':>14}"
        f"{'Phase (deg)':>15}"
    )

    for name, force in zip(
        ["Drag", "Inertia", "Total"], force_arrays
    ):
        i_max = np.argmax(force)
        i_min = np.argmin(force)

        phase_max = np.degrees(omega * t_list[i_max]) % 360
        phase_min = np.degrees(omega * t_list[i_min]) % 360

        print(
            f"{name:<12}"
            f"{force[i_max] / 1000:>14.2f}"
            f"{phase_max:>15.2f}"
            f"{force[i_min] / 1000:>14.2f}"
            f"{phase_min:>15.2f}"
        )


def print_total_contributions(title, drag, inertia):
    total = drag + inertia
    print("\n" + title)

    for label, idx in [
        ("Maximum total", np.argmax(total)),
        ("Minimum total", np.argmin(total))
    ]:
        phase = np.degrees(omega * t_list[idx]) % 360

        print(f"\n{label} — signed values at the SAME instant")
        print(f"Drag:    {drag[idx] / 1000:+.2f} kN")
        print(f"Inertia: {inertia[idx] / 1000:+.2f} kN")
        print(f"Total:   {total[idx] / 1000:+.2f} kN")
        print(f"Time:    {t_list[idx]:.3f} s")
        print(f"Phase:   {phase:.2f} degrees")


print("\nDIAGONAL BRACE — LINEAR WAVE-CURRENT INTERACTION")
print(f"Brace angle:         {np.degrees(alpha):.2f} degrees")
print(f"Current:             {U_current:.3f} m/s along +x")
print(f"Absolute frequency:  {omega:.6f} rad/s")
print(f"Intrinsic frequency: {sigma:.6f} rad/s")
print(f"Wavenumber:          {k:.6f} 1/m")
print(f"Wavelength:          {2 * np.pi / k:.3f} m")
print("Drag: total horizontal velocity and w included to +6 m.")
print("Inertia: ax to 0 m; az to +6 m.")
print("All reported forces are SIGNED.")
print("Phase means omega*t at the front-left column (x=0, z=0).")

print_results(
    "HORIZONTAL FORCES (+x positive)",
    [Fd_x, Fi_x, Ft_x]
)
print_results(
    "VERTICAL FORCES (+z upward)",
    [Fd_z, Fi_z, Ft_z]
)

print_total_contributions(
    "HORIZONTAL TOTAL PEAK DETAILS", Fd_x, Fi_x
)
print_total_contributions(
    "VERTICAL TOTAL PEAK DETAILS", Fd_z, Fi_z
)


# ---------- PLOTS ----------
fig, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=True)
plot_times = np.append(t_list, T)

plot_settings = [
    (axes[0], [Fd_x, Fi_x, Ft_x], "Horizontal forces (+x)"),
    (axes[1], [Fd_z, Fi_z, Ft_z], "Vertical forces (+z upward)")
]

for ax_plot, force_arrays, title in plot_settings:
    for force, name, color in zip(
        force_arrays,
        ["Drag", "Inertia", "Total"],
        ["tab:blue", "tab:purple", "black"]
    ):
        ax_plot.plot(
            plot_times,
            np.append(force, force[0]) / 1000,
            label=name,
            color=color,
            linewidth=2.2 if name == "Total" else 1.4
        )

    total = force_arrays[2]
    extrema = [np.argmax(total), np.argmin(total)]

    ax_plot.scatter(
        t_list[extrema],
        total[extrema] / 1000,
        color="red",
        zorder=5,
        label="Total extrema"
    )

    ax_plot.axhline(0, color="grey", linewidth=0.7)
    ax_plot.set_title(title)
    ax_plot.set_ylabel("Force (kN)")
    ax_plot.grid(alpha=0.3)
    ax_plot.legend()

axes[1].set_xlabel("Time (s)")
fig.suptitle(
    "Diagonal brace — wave-current interaction, current = 1.6 m/s\n"
    "Drag: complete velocity to crest | Inertia: ax to SWL, az to crest"
)

plt.tight_layout()
plt.show()
