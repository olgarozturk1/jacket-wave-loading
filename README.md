# Jacket wave loading

Python calculations and an interactive jacket structure viewer developed for our fluid-loading coursework.

**[Open the interactive structure viewer](https://jacket-structure-viewer.y655375436.chatgpt.site)**

## Naming convention

- `along_x` or `along_y` identifies the beam's orientation.
- `xz` or `yz` identifies the brace's plane.
- `force_x`, `force_y` or `force_z` identifies the force component.
- `from_w` means the contribution caused by vertical wave velocity; `from_u` means horizontal wave velocity.

Every main script includes drag, inertia and their simultaneous total. For the red brace, `forces_x_and_z` means both coordinate directions.

## The eight calculation scripts

The names identify the member orientation, output force direction and velocity contribution. The model assumptions are described below.

| Script | Member | Force output | Velocity dependence |
| --- | --- | --- | --- |
| [beam_along_y_force_x_from_u.py](beam_along_y_force_x_from_u.py) | Lower beam along y | Horizontal x | u and its horizontal acceleration; separate-component drag. |
| [beam_along_x_force_z_from_w.py](beam_along_x_force_z_from_w.py) | Lower beam along x | Vertical z | w and its acceleration; phase varies along x. |
| [beam_along_y_force_z_from_w.py](beam_along_y_force_z_from_w.py) | Lower beam along y | Vertical z | w and its acceleration; constant x. |
| [red_brace_xz_forces_x_and_z.py](red_brace_xz_forces_x_and_z.py) | Brace rising in +x, in the x-z plane | Horizontal x and vertical z | Normal velocity and acceleration projected from u and w. |
| [yellow_brace_yz_force_y_from_w.py](yellow_brace_yz_force_y_from_w.py) | Front yellow brace in the y-z plane | Horizontal y | w-only normal force, resolved into y. |
| [yellow_brace_yz_force_z_from_w.py](yellow_brace_yz_force_z_from_w.py) | Front yellow brace in the y-z plane | Vertical z | w-only normal force, resolved into z. |
| [yellow_brace_yz_force_x_from_u.py](yellow_brace_yz_force_x_from_u.py) | Yellow brace in the y-z plane | Horizontal x | u and its horizontal acceleration; separate-component drag. |
| [vertical_column_force_x_from_u.py](vertical_column_force_x_from_u.py) | Vertical column along z | Horizontal x | u and its acceleration. |

## Run a calculation

Install the required libraries once:

```sh
python -m pip install -r requirements.txt
```

Open a script in Spyder and press Run, or run it from this folder:

```sh
python yellow_brace_yz_force_y_from_w.py
python yellow_brace_yz_force_z_from_w.py
```

Results appear in the console and plots open in a window. Functions return N; tables and plots generally use kN. Total force means drag plus inertia at the same instant, followed by a search for the maximum.

`index.html` is the standalone viewer: download it and open it in a browser. Its sliders illustrate velocity and force projections; they do not calculate wave-load maxima.

## Geometry, coordinates and model status

The viewer follows a 30 m by 30 m plan, 40 m water depth, brace elevations -35 m to +7 m, and a 42 m rise over a 30 m span. The origin is on the front-left column at SWL. +x points towards the back-left column and is the wave direction; +y points towards the front-right column; +z points upward.

The front yellow brace is at x = 0 and falls in +y, from `(0, 0, 7)` to `(0, 30, -35)`. For this geometry, the normal contribution from upward w points towards +y and +z.

### Yellow brace: separate the force caused by w into y and z

`yellow_brace_yz_force_y_from_w.py` reports the horizontal y component; `yellow_brace_yz_force_z_from_w.py` reports the vertical z component. Both drag and inertia are first calculated in the normal direction using the w-only scalar model, then resolved geometrically. Applying the lecture's separate-component treatment to the inclined brace is a modelling assumption.

With alpha = atan(42/30) = 54.4623 degrees, the positive unit normal for the front brace is `(0, sin(alpha), cos(alpha))`. For each of drag, inertia and their simultaneous total:

- `Fy_from_w = Fn_from_w * sin(alpha)` (horizontal, along y).
- `Fz_from_w = Fn_from_w * cos(alpha)` (vertical, along z).
- `Fx_from_w = 0` for this normal-force contribution.

The positive y and z peaks of each force term occur together because they are positive constant multiples of the same signed normal force. Drag, inertia and total have different peak phases. Using the saved fixed-crest parameters at x = 0:

| Force caused by w | Maximum +y (kN) | Maximum +z (kN) | Local phase (degrees) |
| --- | ---: | ---: | ---: |
| Drag | 52.04 | 37.17 | 0.00 |
| Inertia | 20.29 | 14.49 | 270.00 |
| Total at the same instant | 54.02 | 38.58 | 348.76 |

The separate scripts plot y and z individually and print the simultaneous drag/inertia values at maximum total. Their signs match the falling-y front brace in the viewer.

**Phase reference:** the red x-z script still sets x=0 at the brace SWL crossing, so its reference phase differs from the viewer's origin.

### Separate-component drag assumption

The Morisons lecture (printed page 83) gives separate horizontal-cylinder drag terms proportional to `u*abs(u)` and `w*abs(w)`. The horizontal y-beam script uses that horizontal expression directly. For the yellow y-z brace, the calculation applies this separate-component treatment in the two directions normal to its axis:

- The x-force script integrates `0.5*rho*CD*D*u*abs(u)` along the member. The whole horizontal velocity u is normal to this brace; no extra angle projection is needed for the x force.
- The two w-force scripts use `wn = w*cos(alpha)` and calculate `wn*abs(wn)` before resolving the force into y and z.
- This deliberately neglects coupling between u and w in the drag magnitude. It is a stated approximation, not a consequence of the zero x projection of w alone.
- The red x-z brace has both u and w projected into the same normal direction; its signed normal velocity therefore still includes both.

The alternative vector-drag law would use `u*sqrt(u**2+(w*cos(alpha))**2)` for the yellow brace's x drag. In that law w affects the normal-speed magnitude even though its normal projection has no x component. That law is not used in the separate-component scripts here. Inertia is linear; the yellow brace's x inertia depends only on horizontal acceleration in either treatment.

The viewer illustrates projections rather than a drag calculation. Its full normal-speed display can include both u and w without implying that its separate force arrows implement the coupled drag law.

The brace scripts generally use the requested fixed integration limits z = -35 to +6 m for both drag and inertia throughout the cycle. This is a fixed-crest approximation, not an instantaneous wetted-length model. Above-SWL evaluations extend the linear-wave kinematics.

Most scripts use rho = 1025 kg/m³, Cd = 1.2, Cm = 0.8, a = 6 m, k = 0.043 m⁻¹ and omega = 0.628 rad/s. The saved scripts themselves are the source of their parameter values.

The viewer demonstrates geometry and illustrative velocity/force components. Its sliders do not calculate wave-load maxima.

## Working together

Create a branch for a change, explain which assumptions or parameters changed, and compare results using the same coordinates, surface treatment and force component before merging.
