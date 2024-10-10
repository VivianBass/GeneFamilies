import csv

# Archivos de entrada y salida
families_file = 'orthofinder/families.tsv'      # Archivo de familias
orthologs_file = 'orthofinder/orthologs_dmel_dgri.csv'     # Archivo de ortólogos
output_file = 'orthofinder/paralogs_dmel_dgri.tsv'         # Archivo de salida para paralogos

# Leer los ortólogos en un conjunto para una búsqueda rápida
orthologs_set = set()
with open(orthologs_file, 'r', newline='', encoding='utf-8') as f:
    for line in f:
        # Dividir cada línea y añadir los genes a un conjunto
        genes = line.strip().split(', ')
        orthologs_set.update(genes)

# Lista para almacenar los paralogos
paralogs = []

# Leer el archivo de familias
with open(families_file, 'r', newline='', encoding='utf-8') as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
        family_name = row[0]  # Obtener el nombre de la familia
        genes = row[1:]       # Obtener los genes de la familia
        
        for gene in genes:
            # Comprobar si el gen no está en los ortólogos
            if gene not in orthologs_set:
                paralogs.append([f'paralogs_{family_name}', gene])

# Escribir el archivo de salida con los paralogos
with open(output_file, 'w', newline='', encoding='utf-8') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['Family', 'Gene'])  # Escribir la cabecera
    writer.writerows(paralogs)  # Escribir los paralogos

print(f"Archivo '{output_file}' generado con éxito.")
