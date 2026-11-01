import pandas as pd

# Caminhos dos arquivos
input_csv = "/home/esther/Downloads/bioinfo/clinvar_vcv_completo.csv"
output_csv = "/home/esther/Downloads/bioinfo/clinvar_cyp2d6.csv"

print("Carregando o arquivo CSV completo...")
# low_memory=False otimiza a leitura de colunas mistas
df = pd.read_csv(input_csv, low_memory=False)

print("Filtrando registros para o gene CYP2D6...")
# Filtra onde a coluna 'Gene' contém 'CYP2D6' (ignorando maiúsculas/minúsculas e variações)
df_cyp2d6 = df[df['Gene'].str.contains('CYP2D6', na=False, case=False)]

# Salva o resultado filtrado em um novo arquivo menor
df_cyp2d6.to_csv(output_csv, index=False)

print(f"\nFiltragem concluída com sucesso!")
print(f"Total de registros encontrados para o CYP2D6: {len(df_cyp2d6):,}")
print(f"Arquivo salvo em: {output_csv}")

# Exibe uma prévia rápida das primeiras linhas encontradas
print("\n--- Amostra de variantes do CYP2D6 ---")
print(df_cyp2d6[['VCV_Accession', 'Variation_Name', 'Variation_Type', 'Clinical_Significance']].head(10))