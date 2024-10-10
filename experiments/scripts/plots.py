import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

ortologos = pd.read_csv("orthofinder/orthologs_dmel_dgri.tsv", sep='\t')
paralogos = pd.read_csv("orthofinder/paralogs_dmel_dgri.tsv", sep='\t')

conteo_ortologos = ortologos.groupby('Orthogroup').size().reset_index(name='count_ortologos')

conteo_paralogos = paralogos.groupby('Family').size().reset_index(name='count_paralogos')
conteo_paralogos = conteo_paralogos.rename(columns={'Family': 'Orthogroup'})

conteo_total = pd.merge(conteo_ortologos, conteo_paralogos, on='Orthogroup', how='outer').fillna(0)
conteo_total['count_ortologos'] = conteo_total['count_ortologos'].astype(int)
conteo_total['count_paralogos'] = conteo_total['count_paralogos'].astype(int)

total_orths=conteo_total['count_ortologos'].sum()
total_par=conteo_total['count_paralogos'].sum()

grupos_con_mas_de_uno_ortologos = (conteo_total['count_ortologos'] > 1).sum()
grupos_con_mas_de_uno_paralogos = (conteo_total['count_paralogos'] > 1).sum()

prop_orths_par=(grupos_con_mas_de_uno_ortologos/grupos_con_mas_de_uno_paralogos)*100

prop_par_orths=(grupos_con_mas_de_uno_paralogos/grupos_con_mas_de_uno_ortologos)*100

print(f"Orthologs total: {total_orths}")
print(f"Paralogs total: {total_par}")
print(f"Groups with at least two orthologs: {grupos_con_mas_de_uno_ortologos}")
print(f"Groups with at least two paralogs: {grupos_con_mas_de_uno_paralogos}")

print(f"Ratio orth vs paralog: {prop_orths_par:.2f}%")
print(f"Ratio paralog vs orth: {prop_par_orths:.2f}%")

sns.set(style="whitegrid")

plt.figure(figsize=(10, 6))
sns.scatterplot(data=conteo_total, x='count_ortologos', y='count_paralogos', color='blue', s=100)

plt.title("Orthologs vs paralogs per group")
plt.xlabel("Orthologs")
plt.ylabel("Paralogs")

plt.show()



