import psycopg2
from googletrans import Translator
from time import sleep

# Configurações de conexão com o seu banco PostgreSQL
DB_CONFIG = {
    "dbname": "ciera_pgx_db",
    "user": "pregenescan",
    "password": "2Tqm9XnOtk8QMWI2oljq7fnp3jggIvLZ",
    "host": "dpg-db9m7eom7kps73dtm8t0-a.ohio-postgres.render.com",
    "port": "5432",
    "sslmode": "require"
}

NOME_TABELA = "anotacoes_genotipos"

def traduzir_recomendacoes_genotipos():
    conexao = psycopg2.connect(**DB_CONFIG)
    cursor = conexao.cursor()
    translator = Translator()

    try:
        # Busca apenas os registros que ainda não possuem tradução em português
        query = f"""
            SELECT id, recomendacao_clinica_en 
            FROM {NOME_TABELA} 
            WHERE recomendacao_clinica_en IS NOT NULL 
              AND (recomendacao_clinica_pt IS NULL OR recomendacao_clinica_pt = '');
        """
        cursor.execute(query)
        registros = cursor.fetchall()

        total = len(registros)
        print(f"Total de registros pendentes para tradução: {total}")

        if total == 0:
            print("Todos os registros já estão traduzidos!")
            return

        for index, (id_registro, texto_en) in enumerate(registros, start=1):
            if not texto_en or len(texto_en.strip()) < 3:
                continue

            sucesso = False
            tentativas = 0

            while not sucesso and tentativas < 3:
                try:
                    # Realiza a tradução individual com o googletrans
                    resultado = translator.translate(texto_en, src='en', dest='pt')
                    texto_pt = resultado.text

                    if texto_pt:
                        # Atualiza no banco
                        update_query = f"""
                            UPDATE {NOME_TABELA} 
                            SET recomendacao_clinica_pt = %s 
                            WHERE id = %s;
                        """
                        cursor.execute(update_query, (texto_pt, id_registro))
                        conexao.commit()

                        print(f"[{index}/{total}] ID {id_registro} traduzido com sucesso.")
                        sucesso = True
                    else:
                        tentativas += 1
                        sleep(2)

                    # Pausa rápida e segura entre cada registro para evitar bloqueio
                    sleep(0.5)

                except Exception as ex:
                    tentativas += 1
                    print(f"[{index}/{total}] Erro ao traduzir ID {id_registro}: {ex}. Tentativa {tentativas}/3. Aguardando...")
                    sleep(3)

        print(f"\nProcesso de tradução na tabela '{NOME_TABELA}' concluído com sucesso!")

    except Exception as e:
        print(f"Erro geral no banco de dados: {e}")

    finally:
        cursor.close()
        conexao.close()

if __name__ == "__main__":
    traduzir_recomendacoes_genotipos()