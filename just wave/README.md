# Waves only — Q1(a) and Q1(b)

These are the eight wave-only member scripts, moved from the repository root without changing their calculations.

| Script | Member and force component |
| --- | --- |
| [vertical_column_force_x_from_u.py](vertical_column_force_x_from_u.py) | Vertical column: horizontal x force |
| [beam_along_y_force_x_from_u.py](beam_along_y_force_x_from_u.py) | Beam along y: horizontal x force |
| [beam_along_y_force_z_from_w.py](beam_along_y_force_z_from_w.py) | Beam along y: vertical z force |
| [beam_along_x_force_z_from_w.py](beam_along_x_force_z_from_w.py) | Beam along x: vertical z force |
| [red_brace_xz_forces_x_and_z.py](red_brace_xz_forces_x_and_z.py) | Red diagonal in x-z: horizontal x and vertical z forces |
| [yellow_brace_yz_force_x_from_u.py](yellow_brace_yz_force_x_from_u.py) | Yellow diagonal in y-z: horizontal x force from u |
| [yellow_brace_yz_force_y_from_w.py](yellow_brace_yz_force_y_from_w.py) | Yellow diagonal in y-z: horizontal y force from w |
| [yellow_brace_yz_force_z_from_w.py](yellow_brace_yz_force_z_from_w.py) | Yellow diagonal in y-z: vertical z force from w |

Open any script in Spyder and press Run. Each script is standalone. Install the libraries once from this folder with:

```sh
python -m pip install -r requirements.txt
```

For the existing combined wave-only and current-only results summary, run `python part1_results.py` from the repository root.

**Reporting status:** separate drag and inertia extrema are retained in the report. Total-force outputs and their associated phases/times remain provisional working calculations pending assessment of the integration treatment.

See the [main README](../README.md) for the geometry, coordinates and integration assumptions, and the [structure viewer](https://jacket-structure-viewer.y655375436.chatgpt.site/) for the member directions. Wave-and-current calculations are in the separate [wave and current folder](../wave%20and%20current/README.md).
