import os
import json
import pandas as pd

# Directorio que contiene todas las carpetas "cluster_X"
base_dir = "data/phylogeny/"

# Obtener todas las carpetas que comienzan con "cluster_"
carpetas = [d for d in os.listdir(base_dir) if d.startswith("cluster_")]
carpetas.sort()


column_names = ['alpha', 'beta', 'beta-alpha', 'Prob[alpha>beta]', 'Prob[alpha<beta]', 'BayesFactor[alpha<beta]', 'Column_8', 'Column_9']

# Iterar sobre cada carpeta
for carpeta in carpetas:
    print(carpeta)
    # Construir la ruta completa al archivo JSON
    json_path = os.path.join(base_dir, carpeta, f"{carpeta}_CDS_MSA.fasta.FUBAR.json")
    
    # Comprobar si el archivo JSON existe
    if os.path.exists(json_path):
        print(f"Leyendo {json_path}...")
        
        # Leer el archivo JSON
        with open(json_path, 'r') as json_file:
            data = json.load(json_file)
        
        # Extraer la sección que necesitas (por ejemplo, MLE -> content)
        if 'MLE' in data and 'content' in data['MLE']:
            content = data['MLE']['content']['0']
            
            # Convertir a DataFrame
            df = pd.DataFrame(content, columns=column_names)
            
            # Guardar el DataFrame como CSV en la misma carpeta
            csv_path = os.path.join(base_dir, carpeta, f"{carpeta}_ml_tree_no_node_labels.newick.fubar.csv")
            df.to_csv(csv_path, header=False)
            print(f"Archivo CSV guardado en {csv_path}")
        else:
            print(f"Advertencia: El archivo {json_path} no contiene la estructura esperada")
    else:
        print(f"Advertencia: No se encontró el archivo {json_path}")
