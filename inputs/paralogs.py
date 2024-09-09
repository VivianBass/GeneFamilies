import pandas as pd

def txt_to_dataframe(file_path):
    # Leer el archivo .txt
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            # Eliminar espacios en blanco al inicio y final de la línea
            line = line.strip()
            # Separar la línea en "family" y "gene"
            family, gene = line.split()
            # Añadir los datos a la lista
            data.append([family, gene])
    
    # Crear el DataFrame
    df = pd.DataFrame(data, columns=['Family', 'Gene'])
    
    return df

# Ruta al archivo .txt
file_path = 'paralogs.txt'

# Convertir el archivo en un DataFrame
df = txt_to_dataframe(file_path)

print("unique: ", df["Family"].nunique())
# Mostrar el DataFrame
# print(df)

# actual_family=df.iloc[0]["Family"]
# i=0

# for index,row in df.iterrows():
#     current_family=row['Family']
#     if current_family==actual_family:
#         row['Family']="cluster_"+str(i)
#     else:
#         row['Family']="cluster_"+str(i+1)
#         i=i+1
#         actual_family=current_family

# # print(df)
# df.to_csv('paralogs_modify.txt', sep='\t', index=False)
