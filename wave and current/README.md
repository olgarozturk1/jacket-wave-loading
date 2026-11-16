# Q1(d): Wave and current

The eight checked wave-and-current scripts, with clear filenames. Their calculation contents are unchanged from the versions used for the separate drag and inertia results on 15 November 2026. The corrected red-diagonal origin is included.

**Result status:** report the separate drag and inertia extrema. The scripts still calculate and plot simultaneous totals as working calculations, but those total-force maxima and their phases/times remain provisional. They have been left blank in the Word and PDF results tables pending assessment of a consistent wetting/integration treatment.

## Files

| File | Member and force component |
| --- | --- |
| [vertical_column_force_x_from_u.py](vertical_column_force_x_from_u.py) | Vertical column: horizontal x force |
| [beam_along_y_force_x_from_u.py](beam_along_y_force_x_from_u.py) | Beam along y: horizontal x force |
| [beam_along_y_force_z_from_w.py](beam_along_y_force_z_from_w.py) | Beam along y: vertical z force |
| [beam_along_x_force_z_from_w.py](beam_along_x_force_z_from_w.py) | Beam along x: vertical z force |
| [red_diagonal_xz_forces_x_and_z.py](red_diagonal_xz_forces_x_and_z.py) | Red diagonal in x-z: horizontal x and vertical z forces |
| [yellow_diagonal_yz_force_x_from_u.py](yellow_diagonal_yz_force_x_from_u.py) | Yellow diagonal in y-z: horizontal x force from horizontal velocity |
| [yellow_diagonal_yz_force_y_from_w.py](yellow_diagonal_yz_force_y_from_w.py) | Yellow diagonal in y-z: horizontal y force from vertical wave velocity |
| [yellow_diagonal_yz_force_z_from_w.py](yellow_diagonal_yz_force_z_from_w.py) | Yellow diagonal in y-z: vertical z force from vertical wave velocity |

## Run in Spyder

Open any script and press Run. Each is standalone and requires NumPy and Matplotlib. To install the libraries from this folder:

```sh
python -m pip install -r requirements.txt
```

Forces in console tables and plots are in kN. Negative values mean force in the negative coordinate direction.

## Agreed model

- Linear waves with a uniform following current of 1.6 m/s along +x. Local amplitude a = 6 m, depth d = 40 m and absolute frequency omega = 0.628 rad/s are held fixed.
- Solve `(omega - k*Uc)**2 = g*k*tanh(k*d)` on the positive-intrinsic-frequency branch. This gives k = 0.03672575918 per metre and sigma = 0.56923878531 rad/s.
- Wave velocity amplitudes use sigma. Inertia uses linearised material acceleration, including advection by the uniform current. This is the agreed loading-model assumption.
- Signed drag retains `Vn*abs(Vn)`. Beams and the yellow diagonal use separate velocity-component drag; the red diagonal combines u and w in its in-plane normal velocity.
- Diagonal rise/span = 42/30. Yellow length factor is the agreed rounded 1.23; red uses the exact angle factor.

## Coordinates and integration limits

The origin is at the front-left column at SWL; +x follows the waves and current, +y points to the front-right column, and +z is upward. The front yellow diagonal has x = 0. The red diagonal rises from (x,z) = (0,-35) to (30,+7), so `x = (30/42)*(z+35)`.

[Open the structure viewer](https://jacket-structure-viewer.y655375436.chatgpt.site/) for the geometry and component directions.

| Member | Drag limits | Inertia limits |
| --- | --- | --- |
| Column x | z = -40 to +6 m | Horizontal acceleration: -40 to 0 m |
| Yellow diagonal x | z = -35 to +6 m | Horizontal acceleration: -35 to 0 m |
| Yellow diagonal y and z from w | z = -35 to 0 m | Vertical acceleration: -35 to +6 m |
| Red diagonal x and z | Combined normal velocity: -35 to +6 m | Horizontal acceleration: -35 to 0 m; vertical acceleration: -35 to +6 m |
| Beam along y | Full 30 m at z = -35 m | Same member |
| Beam along x | x = 0 to 30 m at z = -35 m | Same member |

These are the preserved fixed, component-specific limits, not an instantaneous moving wetted boundary. Linear kinematics are continued directly above SWL where specified.

Wave-current reference: [Thomas (1981), Wave-current interactions: an experimental and numerical study. Part 1. Linear waves](https://doi.org/10.1017/S0022112081000839).
