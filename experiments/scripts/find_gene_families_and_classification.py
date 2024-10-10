import pandas as pd

experiment="RPKM_paper_diet_m"
fly_or_finder="orthofinder"

families_data = []
with open(fly_or_finder+'/families.tsv', 'r') as file:
    for line in file:
        row = line.strip().split('\t')
        families_data.append(row)

families_df = pd.DataFrame(families_data)
orthologs_df = pd.read_csv(fly_or_finder+'/orthologs.csv')

rpkm_df = pd.read_csv(experiment+'/RPKM.tsv', sep='\t')

ortholog_genes = set(orthologs_df['dmel'])
rpkm_genes = set(rpkm_df['id'])

families_results = {}
orthologs_results = {}
paralogs_results = []

for index, row in families_df.iterrows():
    family_id = row[0]
    genes = row[1:].dropna().tolist()  # List of genes in this family, skiping Nans

    # Filter genes in this family that exists in RPKM 
    genes_in_rpkm = [gene for gene in genes if gene in rpkm_genes]
    
    if genes_in_rpkm:
        
        families_results[family_id] = genes_in_rpkm
        
        
        family_orthologs = [gene for gene in genes_in_rpkm if gene in ortholog_genes]
        family_paralogs = [gene for gene in genes_in_rpkm if gene not in family_orthologs]

        
        if family_orthologs:
            orthologs_results[f'ortholog_{family_id}'] = ','.join(family_orthologs)
        
        
        for gene in family_paralogs:
            paralogs_results.append([f"paralog_{family_id}", gene])

families_output_df = pd.DataFrame.from_dict(families_results, orient='index')
families_output_df.to_csv(experiment+'/families.tsv', sep='\t', header=False)




orthologs_output_df = pd.DataFrame(list(orthologs_results.items()), columns=['family_id', 'dmel'])
orthologs_output_df = orthologs_output_df.drop('family_id', axis=1)
orthologs_output_df.to_csv(experiment+'/orthologs.csv', sep='\t', index=False)

paralogs_output_df = pd.DataFrame(paralogs_results, columns=['Family', 'Gene'])
paralogs_output_df.to_csv(experiment+'/paralogs.csv', sep='\t', index=False)
