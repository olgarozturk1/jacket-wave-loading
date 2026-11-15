"""Run the existing Part 1(a,b) scripts and print one results table.

Also report Part 1(c). Forces are calculated by the member scripts;
this file does not substitute reference values or change their models.
Plot windows are suppressed for this summary run. Run individual member
scripts to display their plots. Requires NumPy and Matplotlib.
"""
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
import runpy

import numpy as np
import matplotlib.pyplot as plt

from part1c_current_only import print_current_results


MEMBERS = [
    ("vertical_column_force_x_from_u.py", "Column Fx"),
    ("beam_along_y_force_x_from_u.py", "Beam along y: Fx"),
    ("yellow_brace_yz_force_x_from_u.py", "Yellow yz: Fx"),
    ("red_brace_xz_forces_x_and_z.py", "Red xz"),
    ("beam_along_y_force_z_from_w.py", "Beam along y: Fz"),
    ("beam_along_x_force_z_from_w.py", "Beam along x: Fz"),
    ("yellow_brace_yz_force_y_from_w.py", "Yellow yz: Fy"),
    ("yellow_brace_yz_force_z_from_w.py", "Yellow yz: Fz"),
]


def collect_wave_results():
    root = Path(__file__).resolve().parent
    rows = []
    for filename, label in MEMBERS:
        previous_show = plt.show
        existing_figures = set(plt.get_fignums())
        try:
            plt.show = lambda *args, **kwargs: None
            with redirect_stdout(StringIO()):
                data = runpy.run_path(str(root / filename), run_name="__main__")
        finally:
            plt.show = previous_show
            for number in set(plt.get_fignums()) - existing_figures:
                plt.close(number)

        if filename in ("yellow_brace_yz_force_y_from_w.py",
                        "yellow_brace_yz_force_z_from_w.py"):
            time = np.linspace(0, 2*np.pi/data["omega"], 100001, endpoint=False)
            fd, fi, total = data["forces"](data["x_brace"], time)
            components = [(label, fd, fi, total, time)]
        elif filename.startswith("red_brace"):
            components = [(label + ": F" + axis, data["Fd_" + axis],
                           data["Fi_" + axis], data["Ft_" + axis], data["t_list"])
                          for axis in ("x", "z")]
        elif "force_z" in filename:
            components = [(label, data["Fd_z"], data["Fi_z"],
                           data["Ftotal_z"], data["time"])]
        else:
            components = [(label, data["Fd"], data["Fi"],
                           data["Ftotal"], data["time"])]

        for name, fd, fi, total, time in components:
            if not np.all(np.isfinite(total)) or not np.allclose(total, fd + fi):
                raise ValueError(f"Invalid simultaneous total in {filename}")
            idx = int(np.argmax(total))
            rows.append(dict(member=name, file=filename,
                             drag_kn=float(np.max(fd)/1000),
                             inertia_kn=float(np.max(fi)/1000),
                             total_kn=float(total[idx]/1000),
                             phase_deg=float(np.degrees(data["omega"]*time[idx]) % 360),
                             time_s=float(time[idx])))
    return rows


def main():
    rows = collect_wave_results()
    print("\nPART 1(a,b): WAVES ONLY — MAXIMUM POSITIVE COMPONENTS")
    print("All forces in kN. Drag and inertia columns are separate maxima.")
    print("Total is maximised AFTER adding signed forces at each instant.")
    print(f"{'Member / component':<24}{'Drag':>11}{'Inertia':>11}"
          f"{'Total':>11}{'Phase (deg)':>14}{'Time (s)':>12}")
    for row in rows:
        print(f"{row['member']:<24}{row['drag_kn']:>11.2f}"
              f"{row['inertia_kn']:>11.2f}{row['total_kn']:>11.2f}"
              f"{row['phase_deg']:>14.2f}{row['time_s']:>12.3f}")
    print("\n+z is upward; front yellow brace: +y and +z normal components.")
    print("Phase is omega*t at each script's x=0 reference.")
    print("Red brace uses its SWL crossing as x=0; see README for coordinates.")
    print("KC numbers are not calculated by these force scripts; see PART1_CHECK.md.")
    print_current_results()
    return rows


if __name__ == "__main__":
    main()
