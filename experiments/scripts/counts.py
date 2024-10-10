import csv

experiment="RPKM_flybase"

# Define los nombres de los archivos de entrada y salida
input_file = experiment+'/families.tsv'
output_file = experiment+'/counts.txt'

# Crear una lista para almacenar los resultados
results = []

# Leer el archivo de entrada
with open(input_file, 'r') as infile:
    reader = csv.reader(infile, delimiter='\t')
    
    # Iterar sobre las filas del archivo
    for row in reader:
        # El primer elemento es el id_group, y el resto son los genes
        id_group = row[0]
        genes = row[1:]
        # print(genes)
        non_empty_genes = [gene for gene in genes if gene]

        # Contar el número de genes
        num_genes = len(non_empty_genes)
        
        # Añadir el resultado a la lista
        results.append([id_group, num_genes])

# Escribir los resultados en el archivo de salida
with open(output_file, 'w', newline='') as outfile:
    writer = csv.writer(outfile, delimiter='\t')
    
    # Escribir el encabezado
    writer.writerow(['id', 'dmel'])
    
    # Escribir los datos
    writer.writerows(results)

print(f'Los conteos se han guardado en {output_file}')
