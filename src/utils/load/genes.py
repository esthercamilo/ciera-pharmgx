import pandas as pd
from sqlalchemy import create_engine

# --- Configurações de Conexão ---
ORIGIN_URL = "postgresql://motiva:kRbQWJ1dmscKoNvIcpdjSlqHilEiqBQl@dpg-dade3agae00c739rl2eg-a.oregon-postgres.render.com/flexreportdb_mmnq"
DEST_URL = "postgresql://postgres:motiva@localhost:5432/ciera_pgx_db"

# Lista de genes desejados
gene_list = ["ABCB1","ABCG2","ACE","ADD1","ADRA2A","ADRB2","AGT","AGTR1","ANKK1","APOA5","APOE","BDKRB1","BDNF","CES1","CNR1","COMT","CYP1A2","CYP2B6","CYP2C9","CYP2C19","CYP2D6","CYP3A4","CYP3A5","CYP4F2","DPYD","DRD1","DRD2","EPHX1","FAAH","F5","FKBP5","G6PD","GNB3","GRIK4","GSK3B","GSTP1","YEATS4","KIF6","HLA-A","HLA-B","HMGCR","HTR1A","HTR2A","HTR2C","MC4R","MTHFR","MTRR","NEDD4L","NOS1AP","OPRD1","OPRM1","PLCD3","RARG","SLC12A3","SLC28A3","SLC6A2","SLC6A4","SLC7A5","SLCO1B1","TPMT","UGT1A1","UGT1A4","UGT1A6","UGT2B15","VKORC1","XPC","rs2952768","CACNA1S","CFTR","MT-RNR1","NAT2","NUDT15","RYR1"]


def copiar_tabela_gene():
    try:
        print("1. Conectando ao banco de origem...")
        engine_origin = create_engine(ORIGIN_URL)

        print("2. Lendo os dados da tabela 'genomics_gene'...")
        # Seleciona id, name, chromosome e pharmgkb_id
        query = "SELECT id, name, chromosome, pharmgkb_id FROM genomics_gene"
        df = pd.read_sql(query, con=engine_origin)

        print(f"-> Total de {len(df)} registros encontrados na origem.")

        print("3. Filtrando os genes desejados...")
        # Filtra as linhas cujo 'name' está na lista de genes
        df = df[df["name"].isin(gene_list)]
        print(f"-> Total de {len(df)} registros após o filtro.")

        # Traduzindo os nomes das colunas de texto para o português (mantendo id e pharmgkb_id)
        df = df.rename(columns={
            "name": "nome",
            "chromosome": "cromossomo"
        })

        print("4. Conectando ao banco de destino...")
        engine_dest = create_engine(DEST_URL)

        print("5. Criando a tabela 'gene' e inserindo os dados no destino...")
        df.to_sql('gene', con=engine_dest, if_exists='replace', index=False)

        print("-> Processo concluído com sucesso! Tabela 'gene' criada e preenchida.")
    except Exception as e:
        print(f"Ocorreu um erro durante a execução: {e}")


if __name__ == "__main__":
    copiar_tabela_gene()