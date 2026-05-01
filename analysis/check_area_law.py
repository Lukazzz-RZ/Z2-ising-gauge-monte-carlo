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
    input_file = "Comprobaciones/LoopsRed24_beta0.72.txt"

print(f"Usando archivo: {input_file}")

n_vals = []
wloop_vals = []
wloop_sigma = []

if not os.path.exists(input_file):
    raise FileNotFoundError(f"No se encontró el fichero: {input_file}")

with open(input_file, 'r') as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split()
        if len(parts) < 3:
            continue
        try:
            n = int(parts[0])
            val = float(parts[1])
            sigma_val = float(parts[2])
            n_vals.append(n)
            wloop_vals.append(val)
            wloop_sigma.append(sigma_val)
        except ValueError:
            continue

if not n_vals:
    raise ValueError("No se han leído datos del fichero.")

# --------------------------
# Transformar eje X (ley de área → n^2)
# --------------------------
x_vals = np.array(n_vals) ** 2
y_vals = np.array(wloop_vals)
y_errs = np.array(wloop_sigma)

# Opcional: ignorar los primeros puntos si son ruidosos
n_fit = 0
x_vals = x_vals[n_fit:]
y_vals = y_vals[n_fit:]
y_errs = y_errs[n_fit:]

# --------------------------
# Ajuste exponencial con errores: W(n) = A * exp(-B * n^2)
# --------------------------
A0 = y_vals[0] if len(y_vals) > 0 else 1.0
B0 = 0.1

popt, pcov = curve_fit(
    exp_func, x_vals, y_vals, p0=[A0, B0],
    sigma=y_errs, absolute_sigma=True  #se propagan los errores reales
)
A_fit, B_fit = popt
perr = np.sqrt(np.diag(pcov))

print("\nAjuste exponencial (ley de área): W(n) = A * exp(-B * n²)")
print(f"A = {A_fit:.4e} ± {perr[0]:.4e}")
print(f"B = {B_fit:.4e} ± {perr[1]:.4e}")

# --------------------------
# Curva ajustada
# --------------------------
x_fit = np.linspace(min(x_vals), max(x_vals), 200)
y_fit = exp_func(x_fit, *popt)

# --------------------------
# Graficar con barras de error
# --------------------------
plt.figure(figsize=(7, 5))
plt.errorbar(x_vals, y_vals, yerr=y_errs, fmt='o', label='Datos', markersize=6, capsize=3)
plt.plot(x_fit, y_fit, '-', label=f"Ajuste: $A e^{{-B n^2}}$\nA={A_fit:.2e}±{perr[0]:.1e}, B={B_fit:.2e}±{perr[1]:.1e}")
plt.xlabel(r'$n^2$')
plt.ylabel(r'$\langle W(\text{loop}) \rangle$')
plt.title('Loops de Wilson promediados (ley de área con errores)')
plt.grid(True)
plt.yscale("log")
plt.legend()
plt.tight_layout()
plt.show()
