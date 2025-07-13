import numpy as np
import csv

sizes = [16, 32, 64, 128, 256]

for n in sizes:
    # Genera matrice casuale di interi 0-255
    mat = np.random.randint(0, 256, size=(n, n), dtype=np.uint8)
    
    # Nome file
    filename = f"matrice_{n}x{n}.csv"
    
    # Salva in CSV
    with open(filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(mat)
    
    print(f"Generata matrice {n}x{n} in file {filename}")
