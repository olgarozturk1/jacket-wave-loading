"""Part 1(c): steady +x current, with no waves.

Current occupies -40 <= z <= 0 m. Brace rise/span are 42/30 m;
its submerged vertical extent is -35 to 0 m. Inertia is zero because
the prescribed uniform current is steady and the members are fixed.
"""
import numpy as np


def current_results(rho=1025.0, CD=1.2, current_speed=1.6):
    """Return horizontal drag, inertia and total in N for four members."""
    depth = 40.0
    beam_length = 30.0
    rise, span = 42.0, 30.0
    submerged_rise = 35.0
    sin_alpha = rise / np.hypot(rise, span)
    brace_length = submerged_rise / sin_alpha
    q = 0.5 * rho * CD

    # (name, diameter, wetted length, normal current, x projection)
    members = [
        ("Vertical column", 4.0, depth, current_speed, 1.0),
        ("Horizontal beam along y", 2.0, beam_length, current_speed, 1.0),
        ("Yellow brace in y-z", 1.0, brace_length, current_speed, 1.0),
        ("Red brace in x-z", 1.0, brace_length,
         current_speed * sin_alpha, sin_alpha),
    ]
    rows = []
    for member, diameter, length, un, normal_x in members:
        drag_x = q * diameter * length * un * abs(un) * normal_x
        rows.append(dict(member=member, diameter_m=diameter,
                         wetted_length_m=length, drag_n=float(drag_x),
                         inertia_n=0.0, total_n=float(drag_x)))
    return rows


def print_current_results(rows=None):
    if rows is None:
        rows = current_results()
    print("\nPART 1(c): CURRENT ONLY — HORIZONTAL FORCE IN +x")
    print("Uniform current = 1.6 m/s; no waves; all forces in kN.")
    print("Column: 40 m wetted; beam: 30 m; brace: 35/sin(alpha) m.")
    print(f"{'Member':<29}{'Drag':>12}{'Inertia':>12}{'Total':>12}")
    for row in rows:
        print(f"{row['member']:<29}{row['drag_n']/1000:>12.3f}"
              f"{row['inertia_n']/1000:>12.3f}{row['total_n']/1000:>12.3f}")
    print("Steady current has no oscillation phase or finite wave-period KC.")


if __name__ == "__main__":
    print_current_results()
