import json
import matplotlib.pyplot as plt
import numpy as np
import sys

json_path = sys.argv[1]
# Carica i dati dal file JSON
with open(json_path) as f:
    data = json.load(f)

# Estrai dati
labels = [item["label"] for item in data["results"]]
N_values = [int(label.split('x')[0]) for label in labels]
homemade_times = [item["homemade_time"] for item in data["results"]]
opencv_times = [item["opencv_time"] for item in data["results"]]

# Grafico 1
plt.figure(figsize=(8,5))
plt.semilogy(N_values, homemade_times, marker='o', label='Homemade DCT')
plt.semilogy(N_values, opencv_times, marker='s', label='OpenCV DCT')
plt.xlabel('Dimensione N (es. lato matrice NxN)')
plt.ylabel('Tempo (s) [scala semilogaritmica]')
plt.title('Confronto tempi DCT2 homemade vs OpenCV')
plt.legend()
#plt.grid(True, which="both", ls="--", linewidth=0.5)
plt.xticks(N_values)
plt.tight_layout()

# Grafico 2
fig, ax = plt.subplots(figsize=(9,6))
x = np.arange(len(N_values))
width = 0.35
ax.bar(x - width/2, homemade_times, width, label='Homemade DCT')
ax.bar(x + width/2, opencv_times, width, label='OpenCV DCT')
ax.set_yscale('log')
ax.set_xlabel('Dimensione N (NxN)')
ax.set_ylabel('Tempo (s) [scala logaritmica]')
ax.set_title('Confronto tempi DCT2 homemade vs OpenCV (istogramma)')
ax.set_xticks(x)
ax.set_xticklabels(N_values)
ax.legend()
#ax.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.tight_layout()

# Mostra tutte le figure **insieme**
plt.show()
