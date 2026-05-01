import matplotlib.pyplot as plt
import numpy as np
import sys
from tkinter import Tk
from tkinter.filedialog import askopenfilename

# Ocultar ventana principal de Tkinter
Tk().withdraw()

# Selección de archivo
if len(sys.argv) > 1:
    filename = sys.argv[1]
else:
    filename = askopenfilename(title="Selecciona el archivo de resultados",
                               filetypes=[("Archivos de texto", "*.txt"), ("Todos los archivos", "*.*")])
    if not filename:
        print("No se seleccionó ningún archivo. Saliendo...")
        exit()

print(f"Usando archivo: {filename}")

# Cargar datos
try:
    data = np.loadtxt(filename)
except Exception as e:
    print(f"Error al leer el fichero '{filename}': {e}")
    exit()

# Separar columnas
attempts = data[:, 0] / 1e6  # Dividir entre 1e6 para que se vea en millones
dE_avg = data[:, 2]

# Filtrar valores válidos
valid = np.isfinite(dE_avg)
dE_valid = dE_avg[valid]
attempts_valid = attempts[valid]

# Tamaños
label_size = 24
tick_size = 24
title_size = 24

# Crear figura
fig, ax = plt.subplots(figsize=(8, 6))

# Graficar ΔE promedio
ax.plot(attempts_valid, dE_valid, marker='o', linestyle=' ', markersize=3, color='red', linewidth=1.2)

# Ejes grandes
ax.set_xlabel(r"Pasos MC $ (\cdot 10^6) $", fontsize=label_size)
ax.set_ylabel("ΔE promedio por bloque", fontsize=label_size)
ax.tick_params(axis='both', labelsize=tick_size)

# Título
ax.set_title(r"$\beta = 0.72$, configuración desordenada", fontsize=title_size)

# Limites de y
ax.set_ylim(dE_valid.min()*1.1, dE_valid.max()*1.1)

# Grid
ax.grid(True, linestyle='--', alpha=0.6)

plt.tight_layout()
plt.show()
