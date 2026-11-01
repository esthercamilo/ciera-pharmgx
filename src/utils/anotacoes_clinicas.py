import csv
import psycopg2
from psycopg2 import Error

# Configurações do seu banco PostgreSQL
DB_CONFIG = {
    "dbname": "ciera_pgx_db",
    "user": "pregenescan",
    "password": "2Tqm9XnOtk8QMWI2oljq7fnp3jggIvLZ",
    "host": "dpg-db9m7eom7kps73dtm8t0-a.ohio-postgres.render.com",
    "port": "5432",
    "sslmode": "require",
}

# Nome do arquivo TSV
TSV_FILE = "bula_do_remedio_FDA.tsv"


def criar_colunas_se_necessario(cursor):
  """Garante que as novas colunas existam na tabela anotacoes_clinicas."""
  queries = [
      (
          "ALTER TABLE anotacoes_clinicas ADD COLUMN IF NOT EXISTS"
          " tem_droga_alternativa VARCHAR(255);"
      ),
      (
          "ALTER TABLE anotacoes_clinicas ADD COLUMN IF NOT EXISTS nivel_do_teste"
          " VARCHAR(255);"
      ),
      (
          "ALTER TABLE anotacoes_clinicas ADD COLUMN IF NOT EXISTS fonte"
          " VARCHAR(255);"
      ),
  ]
  for q in queries:
    cursor.execute(q)


def processar_tabela_pgx():
  sucessos = 0
  falhas = 0
  nao_encontradas = 0
  linhas_total = 0

  print("Conectando ao banco de dados PostgreSQL...")
  try:
    connection = psycopg2.connect(**DB_CONFIG)
    cursor = connection.cursor()

    # 1. Criar colunas novas se não existirem
    criar_colunas_se_necessario(cursor)
    connection.commit()

    print(f"Iniciando a leitura do arquivo {TSV_FILE}...")

    with open(TSV_FILE, mode="r", encoding="utf-8") as file:
      # O delimitador padrão de TSV é a tabulação (\t)
      reader = csv.DictReader(file, delimiter="\t")

      for row in reader:
        linhas_total += 1

        # Extração de campos do TSV conforme mapeamento
        chemical = row.get("Chemicals", "").strip()
        genes_raw = row.get("Genes", "").strip()
        var_haplo = row.get("Variants/Haplotypes", "").strip()
        has_alt_drug = row.get("Has Alternate Drug", "").strip()
        testing_level = row.get("Testing Level", "").strip()
        source = row.get("Source", "").strip()

        # Separar os genes por vírgula ou ponto e vírgula
        # Usamos replace para garantir que vírgulas também sejam tratadas se houverem
        genes_raw = genes_raw.replace(",", ";")
        lista_genes = [g.strip() for g in genes_raw.split(";") if g.strip()]

        # Se não houver genes especificados na linha, tratamos como string vazia para a query
        if not lista_genes:
          lista_genes = [""]

        # Identificar se é haplótipo ou variante (Haplótipo contém asterisco '*')
        haplotipo_id = var_haplo if "*" in var_haplo else None
        variante_id = var_haplo if "*" not in var_haplo and var_haplo else None

        # Para cada gene listado na linha, atualizamos a tabela correspondente
        for gene in lista_genes:
          try:
            # Query de atualização baseada em droga_id (nome_en) e gene_id (nome)
            # Ajuste os nomes das colunas da sua tabela caso sejam diferentes (ex: nome_en, nome)
            update_query = """
                            UPDATE anotacoes_clinicas
                            SET tem_droga_alternativa = %s,
                                nivel_do_teste = %s,
                                fonte = %s,
                                variante_id = COALESCE(%s, variante_id),
                                haplotipo_id = COALESCE(%s, haplotipo_id)
                            WHERE nome_en = %s AND nome = %s;
                        """

            cursor.execute(
                update_query,
                (
                    has_alt_drug,
                    testing_level,
                    source,
                    variante_id,
                    haplotipo_id,
                    chemical,
                    gene,
                ),
            )

            if cursor.rowcount > 0:
              sucessos += 1
            else:
              nao_encontradas += 1

          except Error as e:
            falhas += 1
            print(
                f"Erro na linha {linhas_total} (Droga: {chemical},"
                f" Gene: {gene}): {e}"
            )

      # Salva todas as alterações
      connection.commit()

    print("\n--- RELATÓRIO DE PROCESSAMENTO ---")
    print(f"Total de linhas lidas no TSV: {linhas_total}")
    print(f"Atualizações bem-sucedidas: {sucessos}")
    print(f"Registros não encontrados no banco: {nao_encontradas}")
    print(f"Erros de execução/falhas: {falhas}")

  except Error as e:
    print(f"Erro de conexão ou execução no banco de dados: {e}")
  finally:
    if "connection" in locals() and connection:
      cursor.close()
      connection.close()
      print("Conexão com o banco encerrada.")


if __name__ == "__main__":
  processar_tabela_pgx()