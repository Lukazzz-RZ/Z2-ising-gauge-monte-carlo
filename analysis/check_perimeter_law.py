import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from scipy.optimize import curve_fit

# --------------------------
# Definir función del ajuste exponencial
# --------------------------
def exp_func(x, A, B):
    return A * np.exp(-B * x)

# --------------------------
# Lectura de archivo
# --------------------------
if len(sys.argv) > 1:
    input_file = sys.argv[1]
else:
    input_file = "Results/Comprobaciones/LoopsRed24_beta0.72.txt"

print(f"Usando archivo: {input_file}")
n_vals = []
wloop_vals = []

if not os.path.exists(input_file):
    raise FileNotFoundError(f"No se encontró el fichero: {input_file}")

with open(input_file, 'r') as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            n = int(parts[0])
            val = float(parts[1])
            n_vals.append(n)
            wloop_vals.append(val)
        except ValueError:
            continue

if not n_vals:
    raise ValueError("No se han leído datos del fichero.")

# --------------------------
# Transformar eje X (ley de área → n)
# --------------------------
x_vals = 4*np.array(n_vals) 
y_vals = np.array(wloop_vals)

# --------------------------
# Ajuste exponencial: W(n) = A * exp(-B * 4n)
# --------------------------
A0 = y_vals[0] if len(y_vals) > 0 else 1.0
B0 = 0.1
popt, pcov = popt, pcov = curve_fit(
    exp_func, x_vals, y_vals, p0=[A0, B0],
    method='trf',
    ftol=1e-14, xtol=1e-14, gtol=1e-14, max_nfev=10000
    )
A_fit, B_fit = popt
perr = np.sqrt(np.diag(pcov))

print("\nAjuste exponencial (ley de perimetro): W(n) = A * exp(-B * 4n)")
print(f"A = {A_fit:.4e} ± {perr[0]:.4e}")
print(f"B = {B_fit:.4e} ± {perr[1]:.4e}")

# --------------------------
# Curva ajustada
# --------------------------
x_fit = np.linspace(min(x_vals), max(x_vals), 200)
y_fit = exp_func(x_fit, *popt)

# --------------------------
# Graficar
# --------------------------
plt.figure(figsize=(7, 5))
plt.plot(x_vals, y_vals, 'o', label='Datos', markersize=8)
plt.plot(x_fit, y_fit, '-', label=f"Ajuste: $A e^{{-B 4n}}$\nA={A_fit:.2e}, B={B_fit:.2e}")
plt.xlabel(r'$4n$')
plt.ylabel(r'$\langle W(\text{loop}) \rangle$')
plt.yscale("log")
plt.title('Loops de Wilson promediados con ajuste exponencial (ley de perimetro)')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
