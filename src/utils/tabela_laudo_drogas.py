import csv
import psycopg2

# Configurações do seu banco PostgreSQL
DB_CONFIG = {
    "dbname": "ciera_pgx_db",
    "user": "pregenescan",
    "password": "2Tqm9XnOtk8QMWI2oljq7fnp3jggIvLZ",
    "host": "dpg-db9m7eom7kps73dtm8t0-a.ohio-postgres.render.com",
    "port": "5432",
    "sslmode": "require"
}


def importar_gene_medicamento():
    # 1. Ler o CSV e criar um set único de medicamentos
    medicamentos_set = set()

    try:
        # utf-8-sig remove automaticamente caracteres invisíveis de BOM que o Excel costuma salvar no CSV
        with open('gene_medicamento.csv', mode='r', encoding='utf-8-sig') as arquivo_csv:
            leitor = csv.DictReader(arquivo_csv)
            for linha in leitor:
                med = linha.get('medicamento', '').strip()
                if med:
                    medicamentos_set.add(med)

        print(f"Total de medicamentos únicos extraídos do CSV: {len(medicamentos_set)}")

    except FileNotFoundError:
        print("Erro: O arquivo 'gene_medicamento.csv' não foi encontrado na mesma pasta do script.")
        return

    # 2. Conectar ao PostgreSQL
    conexao = psycopg2.connect(**DB_CONFIG)
    cursor = conexao.cursor()

    try:
        laudo_id = 1
        inseridos = 0
        nao_encontrados = []

        print("Iniciando inserção na tabela 'laudo_drogas'...")

        for medicamento in medicamentos_set:
            # Busca o ID na tabela 'drogas' procurando tanto em nome_pt quanto em nome_en
            cursor.execute("""
                SELECT id FROM drogas 
                WHERE nome_pt ILIKE %s OR nome_en ILIKE %s
                LIMIT 1;
            """, (medicamento, medicamento))

            resultado = cursor.fetchone()

            if resultado:
                droga_id = resultado[0]

                # Insere na tabela de ligação.
                # O 'ON CONFLICT DO NOTHING' evita erros caso a relação laudo_id + droga_id já exista.
                cursor.execute("""
                    INSERT INTO laudo_drogas (laudo_id, droga_id)
                    VALUES (%s, %s)
                    ON CONFLICT DO NOTHING;
                """, (laudo_id, droga_id))

                inseridos += 1
            else:
                nao_encontrados.append(medicamento)

        conexao.commit()
        print(f"\nProcesso finalizado com sucesso!")
        print(f"- Relações inseridas/verificadas: {inseridos}")

        if nao_encontrados:
            print(f"- Atenção: {len(nao_encontrados)} medicamentos do CSV não foram encontrados na tabela 'drogas':")
            for med in nao_encontrados:
                print(f"   * {med}")

    except Exception as e:
        conexao.rollback()
        print(f"Erro durante a execução no banco de dados: {e}")

    finally:
        cursor.close()
        conexao.close()


if __name__ == "__main__":
    importar_gene_medicamento()