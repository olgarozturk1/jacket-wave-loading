import numpy as np
import matplotlib.pyplot as plt

# ---------- PARAMETERS ----------
rho = 1025.0
Cd = 1.2
Cm = 0.8
D = 1.0

a = 6.0
k = 0.043
omega = 0.628
d = 40.0

z_min = -35.0
z_swl = 0.0
z_crest = 6.0

rise = 42.0
span = 30.0

alpha = np.arctan2(rise, span)
s = np.sin(alpha)
c = np.cos(alpha)
m = span / rise

T = 2 * np.pi / omega

# Brace rises in +x.
# x = 0 at the SWL crossing, so x = m*z.
#
# Drag:
#   Both u and w are included from -35 m to +6 m.
#
# Inertia — previous prescribed cutoffs retained:
#   ax included from -35 m to 0 m.
#   az included from -35 m to +6 m.


# ---------- DEPTH AND TIME GRIDS ----------
z_lower = np.linspace(z_min, z_swl, 4001)
z_upper = np.linspace(z_swl, z_crest, 2001)

t_list = np.linspace(0, T, 8000, endpoint=False)


def integrate(values, z):
    """Trapezoidal integration with respect to z."""
    return np.sum(
        0.5 * (values[:-1] + values[1:]) * np.diff(z)
    )


B = a * omega / np.sinh(k * d)

cosh_lower = np.cosh(k * (z_lower + d))
sinh_lower = np.sinh(k * (z_lower + d))

cosh_upper = np.cosh(k * (z_upper + d))
sinh_upper = np.sinh(k * (z_upper + d))

Fd_x = np.zeros(t_list.size)
Fi_x = np.zeros(t_list.size)
Fd_z = np.zeros(t_list.size)
Fi_z = np.zeros(t_list.size)


# ---------- CALCULATE FORCES ----------
for i, t in enumerate(t_list):

    # ----- BELOW SWL: -35 m TO 0 m -----
    theta_lower = omega * t - k * m * z_lower

    u_lower = B * cosh_lower * np.sin(theta_lower)
    w_lower = B * sinh_lower * np.cos(theta_lower)

    ax_lower = B * omega * cosh_lower * np.cos(theta_lower)
    az_lower = -B * omega * sinh_lower * np.sin(theta_lower)

    # Complete normal velocity
    un_lower = u_lower * s - w_lower * c

    # Both acceleration components included below SWL
    an_lower = ax_lower * s - az_lower * c

    # ----- ABOVE SWL: 0 m TO +6 m -----
    theta_upper = omega * t - k * m * z_upper

    u_upper = B * cosh_upper * np.sin(theta_upper)
    w_upper = B * sinh_upper * np.cos(theta_upper)

    # CORRECTION: retain BOTH u and w above SWL
    un_upper = u_upper * s - w_upper * c

    # Inertia: retain previous cutoff for ax at SWL
    az_upper = -B * omega * sinh_upper * np.sin(theta_upper)
    an_upper = -az_upper * c

    # ----- INTEGRATE NORMAL LOADS -----
    I_drag = (
        integrate(un_lower * np.abs(un_lower), z_lower)
        + integrate(un_upper * np.abs(un_upper), z_upper)
    )

    I_inertia = (
        integrate(an_lower, z_lower)
        + integrate(an_upper, z_upper)
    )

    # Actual member length: dl = dz / sin(alpha)
    Fd_normal = 0.5 * rho * Cd * D * I_drag / s
    Fi_normal = rho * Cm * (np.pi * D**2 / 4) * I_inertia / s

    # Horizontal components: positive in +x
    Fd_x[i] = Fd_normal * s
    Fi_x[i] = Fi_normal * s

    # Vertical components: positive upwards
    Fd_z[i] = -Fd_normal * c
    Fi_z[i] = -Fi_normal * c


# ---------- TOTAL FORCES AT THE SAME TIME ----------
Ft_x = Fd_x + Fi_x
Ft_z = Fd_z + Fi_z


# ---------- PEAK RESULTS ----------
def print_results(title, force_arrays, positive=True):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)

    print(
        f"{'Force':<15}"
        f"{'Magnitude (kN)':>17}"
        f"{'Time (s)':>16}"
        f"{'Phase (deg)':>18}"
    )
    print("-" * 72)

    for name, force in zip(
        ["Drag", "Inertia", "TOTAL"], force_arrays
    ):
        idx = np.argmax(force) if positive else np.argmin(force)

        magnitude = abs(force[idx]) / 1000
        time = t_list[idx]
        phase = np.degrees(omega * time) % 360

        print(
            f"{name:<15}"
            f"{magnitude:>17.2f}"
            f"{time:>16.3f}"
            f"{phase:>18.2f}"
        )


def print_total_contributions(title, drag, inertia):
    total = drag + inertia

    print("\n" + title)

    for label, idx in [
        ("Maximum positive total", np.argmax(total)),
        ("Maximum negative total", np.argmin(total))
    ]:
        phase = np.degrees(omega * t_list[idx]) % 360

        print(f"\n{label} — signed values at the SAME instant")
        print(f"Drag:    {drag[idx] / 1000:+.2f} kN")
        print(f"Inertia: {inertia[idx] / 1000:+.2f} kN")
        print(f"Total:   {total[idx] / 1000:+.2f} kN")
        print(f"Time:    {t_list[idx]:.3f} s")
        print(f"Phase:   {phase:.2f} degrees")


print("\nDIAGONAL BRACE — WAVE FORCE RESULTS")
print(f"Brace angle: {np.degrees(alpha):.2f} degrees")
print(f"Water depth: {d:.1f} m")
print(f"Lower integration limit: {z_min:.1f} m")
print("Drag: both u and w included up to +6 m.")
print("Inertia: ax included to 0 m; az included to +6 m.")
print("All peak-table entries are positive force MAGNITUDES.")

print_results(
    "HORIZONTAL — PEAK FORCE TOWARDS +x",
    [Fd_x, Fi_x, Ft_x],
    positive=True
)

print_results(
    "HORIZONTAL — PEAK FORCE TOWARDS -x",
    [Fd_x, Fi_x, Ft_x],
    positive=False
)

print_results(
    "VERTICAL — PEAK UPWARD FORCE",
    [Fd_z, Fi_z, Ft_z],
    positive=True
)

print_results(
    "VERTICAL — PEAK DOWNWARD FORCE",
    [Fd_z, Fi_z, Ft_z],
    positive=False
)

print_total_contributions(
    "HORIZONTAL TOTAL PEAK DETAILS", Fd_x, Fi_x
)

print_total_contributions(
    "VERTICAL TOTAL PEAK DETAILS", Fd_z, Fi_z
)

print("\nTOTAL is maximised after adding forces at each time.")
print("It is NOT the sum of separate drag and inertia maxima.")
print("Phase means omega*t at the reference position x = 0.")


# ---------- PLOTS ----------
fig, axes = plt.subplots(
    2, 1, figsize=(11, 8), sharex=True
)

plot_times = np.append(t_list, T)

plot_settings = [
    (axes[0], [Fd_x, Fi_x, Ft_x], "Horizontal forces (+x)"),
    (axes[1], [Fd_z, Fi_z, Ft_z], "Vertical forces (+z upwards)")
]

for ax, force_arrays, title in plot_settings:

    for force, name, color in zip(
        force_arrays,
        ["Drag", "Inertia", "Total"],
        ["tab:blue", "tab:purple", "black"]
    ):
        plot_force = np.append(force, force[0]) / 1000

        ax.plot(
            plot_times,
            plot_force,
            label=name,
            color=color,
            linewidth=2.2 if name == "Total" else 1.4
        )

    total = force_arrays[2]

    for j, idx in enumerate([
        np.argmax(total), np.argmin(total)
    ]):
        ax.scatter(
            t_list[idx],
            total[idx] / 1000,
            color="red",
            zorder=5,
            label="Total peaks" if j == 0 else None
        )

    ax.axhline(0, color="grey", linewidth=0.7)
    ax.set_title(title)
    ax.set_ylabel("Force (kN)")
    ax.grid(alpha=0.3)
    ax.legend()

axes[1].set_xlabel("Time (s)")

fig.suptitle(
    "Diagonal brace — drag includes u and w to crest\n"
    "Inertia: ax to SWL; az to crest"
)

plt.tight_layout()
plt.show()
