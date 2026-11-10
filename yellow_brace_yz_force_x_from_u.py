"""Yellow y-z brace: horizontal x forces from u, using separate-component drag.
Lecture-style assumption: calculate u and w normal-force contributions separately.
Fixed crest integration z=-35 to +6 m. Requested rounded dl/dz=1.23.
The brace is at x=0, so reference phase equals local phase everywhere.
"""
import numpy as np
import matplotlib.pyplot as plt

rho, CD, CM = 1025.0, 1.2, 0.8
D, a, d = 1.0, 6.0, 40.0
k, omega = 0.043, 0.628
length_factor = 1.23

z = np.linspace(-35.0, 6.0, 4001)
weights = np.full(z.size, z[1] - z[0])
weights[[0, -1]] *= 0.5
B = a * omega / np.sinh(k*d)
U = B * np.cosh(k*(z+d))

phase = np.linspace(0, 2*np.pi, 8001)
time = phase / omega
Fd = np.zeros_like(phase)
for j, theta in enumerate(phase):
    u = U * np.sin(theta)
    # u is entirely normal to the y-z brace and acts along x.
    Fd[j] = length_factor * 0.5*rho*CD*D * np.sum(weights*u*np.abs(u))

Fi = (length_factor * rho*CM*np.pi*D**2/4 * omega
      * np.sum(weights*U) * np.cos(phase))
Ftotal = Fd + Fi

print('\nYELLOW BRACE: HORIZONTAL x FORCES FROM u')
print('Fixed limits z=-35 to +6 m; dl/dz=1.23; x=0')
print('\nForce       Maximum positive (kN)    Phase (degrees)')
for name, force in [('Drag', Fd), ('Inertia', Fi), ('Total', Ftotal)]:
    j = np.argmax(force)
    print(f'{name:<12}{force[j]/1000:>18.2f}{np.degrees(phase[j]):>20.2f}')
j = np.argmax(Ftotal)
print(f'\nAt total maximum: drag={Fd[j]/1000:.2f} kN, '
      f'inertia={Fi[j]/1000:.2f} kN, time={time[j]:.3f} s')

plt.figure(figsize=(10, 5))
for name, force in [('Drag', Fd), ('Inertia', Fi), ('Total', Ftotal)]:
    plt.plot(np.degrees(phase), force/1000, label=name)
plt.scatter(np.degrees(phase[j]), Ftotal[j]/1000, color='red', zorder=5)
plt.axhline(0, color='black', linewidth=0.7)
plt.xlabel('Wave phase at x=0 (degrees)')
plt.ylabel('Horizontal force (kN)')
plt.title('Yellow y-z brace: horizontal x forces from u; fixed crest limits')
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()
