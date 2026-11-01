import pandas as pd
from sqlalchemy import create_engine, text

# --- Configurações de Conexão ---
ORIGIN_URL = "postgresql://postgres:motiva@localhost:5432/flexreportdb"
DEST_URL = "postgresql://postgres:motiva@localhost:5432/ciera_pgx_db"


def copiar_tabela_haplotipos():
    try:
        print("1. Conectando aos bancos de dados (origem e destino)...")
        engine_origin = create_engine(ORIGIN_URL)
        engine_dest = create_engine(DEST_URL)

        print("2. Lendo os IDs válidos de genes da tabela 'gene' no destino...")
        # Lê apenas a coluna 'id' da tabela gene existente no destino
        df_genes_destino = pd.read_sql("SELECT id FROM gene", con=engine_dest)
        valid_gene_ids = set(df_genes_destino["id"])
        print(f"-> Encontrados {len(valid_gene_ids)} genes válidos no destino.")

        print("3. Lendo os dados da tabela 'genomics_haplotype' na origem...")
        # Seleciona as colunas necessárias da origem
        query = "SELECT id, name, description, function, link, gene_id FROM genomics_haplotype"
        df_haplotypes = pd.read_sql(query, con=engine_origin)
        print(f"-> Total de {len(df_haplotypes)} registros encontrados na origem.")

        print("4. Filtrando haplótipos para manter apenas os gene_ids existentes na tabela 'gene'...")
        # Filtra mantendo apenas as linhas onde gene_id está no conjunto de IDs válidos
        df_filtered = df_haplotypes[df_haplotypes["gene_id"].isin(valid_gene_ids)].copy()
        print(f"-> Total de {len(df_filtered)} registros após o filtro.")

        print("5. Traduzindo os nomes das colunas para o português...")
        # Mapeia name -> nome, description -> descricao e function -> funcao
        df_filtered = df_filtered.rename(columns={
            "name": "nome",
            "description": "descricao",
            "function": "funcao"
        })

        print("6. Criando a tabela 'haplotipos' e inserindo os dados no destino...")
        # if_exists='replace' substitui a tabela 'haplotipos' caso ela já exista
        df_filtered.to_sql('haplotipos', con=engine_dest, if_exists='replace', index=False)

        print("-> Processo concluído com sucesso! Tabela 'haplotipos' criada e preenchida.")

    except Exception as e:
        print(f"Ocorreu um erro durante a execução: {e}")

def criar_e_popular_haplotipos_variantes():
    try:
        print("1. Conectando aos bancos de dados (origem e destino)...")
        engine_origin = create_engine(ORIGIN_URL)
        engine_dest = create_engine(DEST_URL)

        print("2. Garantindo chave primaria na tabela 'haplotipos'...")
        with engine_dest.begin() as conn:
            # Adiciona a chave primaria em haplotipos caso ainda nao tenha
            conn.execute(text("""
                DO $$ 
                BEGIN 
                    IF NOT EXISTS (
                        SELECT 1 FROM pg_constraint c
                        JOIN pg_namespace n ON n.oid = c.connamespace
                        WHERE c.contype = 'p' AND c.conrelid = 'haplotipos'::regclass
                    ) THEN
                        ALTER TABLE haplotipos ADD PRIMARY KEY (id);
                    END IF;
                END $$;
            """))

            print("3. Criando a tabela 'haplotipos_variantes' no destino...")
            conn.execute(text("""
                CREATE TABLE IF NOT EXISTS haplotipos_variantes (
                    id INT PRIMARY KEY,
                    haplotipo_id INT NOT NULL,
                    variante_id INT NOT NULL,
                    CONSTRAINT fk_haplotipo_variante_haplotipo 
                        FOREIGN KEY (haplotipo_id) 
                        REFERENCES haplotipos(id) 
                        ON DELETE CASCADE,
                    CONSTRAINT fk_haplotipo_variante_variante 
                        FOREIGN KEY (variante_id) 
                        REFERENCES variantes(id) 
                        ON DELETE CASCADE
                );
            """))

        print("4. Lendo os IDs válidos de haplotipos e variantes no destino...")
        df_haplotipos_dest = pd.read_sql("SELECT id FROM haplotipos", con=engine_dest)
        valid_haplotype_ids = set(df_haplotipos_dest["id"])

        df_variantes_dest = pd.read_sql("SELECT id FROM variantes", con=engine_dest)
        valid_variant_ids = set(df_variantes_dest["id"])

        print(f"-> Encontrados {len(valid_haplotype_ids)} haplótipos e {len(valid_variant_ids)} variantes válidos no destino.")

        print("5. Lendo os dados da tabela 'genomics_haplotype_variants' na origem...")
        query = "SELECT id, haplotype_id, variant_id FROM genomics_haplotype_variants"
        df_hv = pd.read_sql(query, con=engine_origin)
        print(f"-> Total de {len(df_hv)} registros encontrados na origem.")

        print("6. Filtrando para manter apenas associações com haplótipos e variantes válidos...")
        df_filtered = df_hv[
            df_hv["haplotype_id"].isin(valid_haplotype_ids) &
            df_hv["variant_id"].isin(valid_variant_ids)
        ].copy()
        print(f"-> Total de {len(df_filtered)} registros após o filtro.")

        print("7. Traduzindo os nomes das colunas para o português...")
        df_filtered = df_filtered.rename(columns={
            "haplotype_id": "haplotipo_id",
            "variant_id": "variante_id"
        })

        print("8. Inserindo os dados na tabela 'haplotipos_variantes' no destino...")
        df_filtered.to_sql('haplotipos_variantes', con=engine_dest, if_exists='append', index=False)

        print("9. Criando os índices para otimização de busca...")
        with engine_dest.begin() as conn:
            conn.execute(text("CREATE INDEX IF NOT EXISTS idx_haplotipos_variantes_haplotipo_id ON haplotipos_variantes(haplotipo_id);"))
            conn.execute(text("CREATE INDEX IF NOT EXISTS idx_haplotipos_variantes_variante_id ON haplotipos_variantes(variante_id);"))

        print("Processo concluído com sucesso! Tabela 'haplotipos_variantes' populada e índices criados.")

    except Exception as e:
        print(f"Ocorreu um erro durante a execução: {e}")


if __name__ == "__main__":
    # copiar_tabela_haplotipos()
    criar_e_popular_haplotipos_variantes()