# Part 1 calculation check

Checked on 15 November 2026 against the group comparison document. The eight member scripts were executed from GitHub commit `c95d6527689bdef10db46eb11126e5ecba366fdf`. The group confirmed that the total-force and phase columns still needed updating. The existing member formulas and prescribed integration limits were retained during this check.

Run `python part1_results.py` to calculate the wave-force table and the current-only results from the repository scripts. No reference force values are substituted into the calculations.

## Parts 1(a) and 1(b)

All forces below are in kN. Each drag or inertia value is its own positive maximum. Total is the maximum of the signed sum at the same instant. Phase is the reference phase at maximum positive total, in degrees.

| Member / component | Drag | Inertia | Updated total | Updated phase | Time (s) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Column Fx | 835.91 | 567.06 | 932.08 | 70.16 | 1.950 |
| Beam along y: Fx | 75.09 | 69.23 | 91.05 | 62.55 | 1.738 |
| Yellow yz: Fx | 249.58 | 40.10 | 251.19 | 85.41 | 2.374 |
| Red xz: Fx | 165.57 | 28.34 | 162.55 | 109.94 | 3.055 |
| Beam along y: Fz | 3.37 | 14.66 | 14.66 | 270.00 | 7.504 |
| Beam along x: Fz | 2.94 | 13.66 | 13.80 | 315.00 | 8.754 |
| Yellow yz: Fy | 28.49 | 20.29 | 32.10 | 339.14 | 9.425 |
| Yellow yz: Fz | 20.35 | 14.49 | 22.93 | 339.14 | 9.425 |
| Red xz: Fz | 118.27 | 20.24 | 116.11 | 289.94 | 8.058 |

The front yellow brace falls in +y. Its positive y and z total peaks occur at the same phase. The red-brace script uses the brace SWL crossing as its x=0 phase reference, whereas the viewer uses the front-left column. Its phase values must not be combined with other members for global loads until that reference is reconciled.

## Comparison with the existing table

| Member / component | Table total | Calculated total | Table phase | Calculated phase |
| --- | ---: | ---: | ---: | ---: |
| Column Fx | 1002.40 | 932.08 | 63.59 | 70.16 |
| Beam along y: Fx | 91.13 | 91.05 | 62.60 | 62.55 |
| Yellow yz: Fx | 252.89 | 251.19 | 82.94 | 85.41 |
| Red xz: Fx | 160.45 | 162.55 | 108.14 | 109.94 |
| Beam along y: Fz | 14.66 | 14.66 | 270.00 | 270.00 |
| Beam along x: Fz | 13.80 | 13.80 | 315.00 | 315.00 |
| Yellow yz: Fy | 54.02 | 32.10 | 168.76 | 339.14 |
| Yellow yz: Fz | 38.58 | 22.93 | 348.76 | 339.14 |
| Red xz: Fz | 114.61 | 116.11 | 288.13 | 289.94 |

The newer cutoffs change both total maxima and their phases. Two original table entries demonstrate why totals needed refreshing: 54.02 kN exceeds 28.49 + 20.29 = 48.78 kN for yellow Fy, and 38.58 kN exceeds 20.35 + 14.49 = 34.84 kN for yellow Fz. A simultaneous total cannot exceed the sum of the separate positive maxima.

### Remaining drag and inertia rounding differences

- Column drag: table 836.8 kN; saved script 835.91 kN (about 0.106%). Using omega = 2*pi/10 instead of 0.628 gives 836.76 kN, which rounds to 836.8. However, that same change also changes column inertia from 567.06 to 567.63 kN. The table currently mixes these precision choices.
- Horizontal y-beam Fx: table 75.16 / 69.3 kN; saved script 75.09 / 69.23 kN. Using omega = 2*pi/10 gives 75.163 / 69.303 kN and explains the table values.
- The other wave drag and inertia entries agree to the displayed two decimal places. No member coefficient or amplitude was tuned to match a rounded table entry. All saved member scripts retain omega = 0.628 rad/s and k = 0.043 1/m.

## Prescribed integration limits actually used

| Member / contribution | Drag limits | Inertia limits |
| --- | --- | --- |
| Column, from u | z = -40 to +6 m | z = -40 to 0 m |
| Yellow yz, x force from u | z = -35 to +6 m | z = -35 to 0 m |
| Yellow yz, y and z forces from w | z = -35 to 0 m | z = -35 to +6 m |
| Red xz, both normal-velocity contributions | z = -35 to +6 m | ax contribution: -35 to 0 m; az contribution: -35 to +6 m |
| Lower beam along y | 30 m length at z = -35 m | Same member |
| Lower beam along x | x = 0 to 30 m, at z = -35 m | Same member |

These are the prescribed component-specific approximations in the saved coursework calculations. They do not describe a single instantaneous wetted boundary. Above-SWL evaluations extend the linear-wave kinematics. This audit checks implementation and consistency with those choices; it does not independently validate the choices.

## Part 1(c) current only

`part1c_current_only.py` implements the four group-approved current-only targets, using rho = 1025 kg/m^3, Cd = 1.2, uniform current = 1.6 m/s, water depth = 40 m, and brace rise/span = 42/30 m. The actual submerged brace length is 35/sin(alpha).

| Member | Comparison table horizontal drag (kN) | Calculation (kN) | Inertia (kN) |
| --- | ---: | ---: | ---: |
| Column | 251.904 | 251.904000 | 0 |
| Beam along y | 94.464 | 94.464000 | 0 |
| Yellow yz brace | 67.717 | 67.717505 | 0 |
| Red xz brace | 36.487 | 36.487784 | 0 |

All four agree within 0.001 kN. Rounding to three decimal places gives 67.718 and 36.488 kN for the braces; the existing table values appear truncated at that precision. Total equals drag for this steady uniform current. No wave phase applies.

## KC coverage and definitions

The eight wave-force scripts do not calculate or print KC. Therefore the KC column is not verified merely by running them. KC requires a stated velocity amplitude and reference elevation. With KC = velocity_amplitude*T/D, the saved wave parameters, and T = 2*pi/omega:

| Reference velocity | KC | Relation to the table |
| --- | ---: | --- |
| Column, horizontal U at SWL, D=4 m | 10.0492 | 10.05 |
| Lower beam, horizontal U at z=-35 m, D=2 m | 7.1361 | 7.14 |
| Lower beam, vertical W at z=-35 m, D=2 m | 1.5110 | 1.51 |
| Yellow yz, horizontal U at SWL, D=1 m | 40.1968 | 40.2 |
| Red xz, horizontal-only projection U*sin(alpha) at SWL | 32.7095 | 32.7 |
| Projection U*cos(alpha) at SWL | 23.3639 | Numerically reproduces 23.4, but uses horizontal U |
| Yellow yz, actual w-only normal amplitude W*cos(alpha) at SWL | 21.9122 | Differs from 23.4 |
| Red xz, full normal amplitude sqrt((U*sin(alpha))^2+(W*cos(alpha))^2) at SWL | 39.3707 | Differs from horizontal-only 32.7 |

The last distinctions need an explicit reference-velocity convention in the report. The confirmed KC entries have not been relabelled or silently replaced. KC is a velocity-based quantity, not a sum of drag and inertia forces.

## Verification

All eight original scripts ran with their plotting code using a non-interactive backend. All nine force components were finite and their total arrays equalled drag plus inertia at every time sample. The combined summary was checked against these independent script outputs. Current-only forces were also checked against the closed-form projected Morison drag integrals.
