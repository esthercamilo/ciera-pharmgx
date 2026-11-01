P = "2Tqm9XnOtk8QMWI2oljq7fnp3jggIvLZ"

import csv
import psycopg2

DB_CONFIG = {
    "dbname": "ciera_pgx_db",
    "user": "pregenescan",
    "password": P,
    "host": "dpg-db9m7eom7kps73dtm8t0-a.ohio-postgres.render.com",
    "port": "5432",
    "sslmode": "require"
}

CSV_INPUT = "/home/esther/GitProjects/ciera-pharmgx/src/utils/load/genes_rs_completo.csv"
BED_OUTPUT = "/home/esther/GitProjects/ciera-pharmgx/src/utils/load/variantes_ampliadas.bed"


def gerar_bed_de_variantes(config, csv_path, output_bed_path):
    conn = None
    cursor = None
    registros_gerados = 0
    total_linhas = 0

    try:
        conn = psycopg2.connect(**config)
        cursor = conn.cursor()

        with open(csv_path, mode='r', encoding='utf-8-sig') as infile, \
                open(output_bed_path, mode='w', encoding='utf-8') as outfile:

            leitor = csv.DictReader(infile, delimiter=',')
            print(f"-> Cabeçalhos lidos: {leitor.fieldnames}")

            for linha in leitor:
                linha_limpa = {k.strip().lower(): v.strip() for k, v in linha.items() if k}

                gene = linha_limpa.get('gene', '')
                rs = linha_limpa.get('rs', '')

                if not gene or not rs:
                    continue

                total_linhas += 1

                query = """
                    SELECT 
                        g.cromossomo AS chrom, 
                        v.start, 
                        v.stop, 
                        v.nome, 
                        g.nome AS gene_nome
                    FROM variantes v
                    JOIN genes g ON v.gene_id = g.id
                    WHERE v.nome ILIKE %s AND v.start IS NOT NULL AND v.stop IS NOT NULL;
                """

                cursor.execute(query, (rs,))
                resultados = cursor.fetchall()

                if not resultados:
                    print(f"⚠️ [Não encontrado ou sem coordenadas] Gene: {gene} | RS: {rs}")
                    continue

                for chrom, v_start, v_stop, nome, gene_nome in resultados:
                    # Regra de ajuste de coordenadas:
                    # Se start == stop, define stop como stop + 1
                    # Se forem diferentes, mantém os valores originais
                    if v_start == v_stop:
                        start_adj = v_start
                        stop_adj = v_stop + 1
                    else:
                        start_adj = v_start
                        stop_adj = v_stop

                    label = f"{nome}_{gene_nome}" if gene_nome else nome

                    chrom_str = str(chrom).strip() if chrom else "chr1"

                    if chrom_str in ['chrM', 'M', 'chrMT']:
                        chrom_str = 'MT'
                    elif not chrom_str.startswith('chr') and chrom_str != 'MT':
                        chrom_str = f"chr{chrom_str}"

                    outfile.write(f"{chrom_str}\t{start_adj}\t{stop_adj}\t{label}\n")
                    registros_gerados += 1

        print(f"\n========================================")
        print(f"Processamento concluído!")
        print(f"Total de linhas lidas no CSV: {total_linhas}")
        print(f"Regiões válidas escritas no BED: {registros_gerados}")
        print(f"Arquivo gerado em: {output_bed_path}")
        print(f"========================================")

    except Exception as e:
        print(f"Erro: {e}")
    finally:
        if cursor:
            cursor.close()
        if conn:
            conn.close()


if __name__ == "__main__":
    gerar_bed_de_variantes(DB_CONFIG, CSV_INPUT, BED_OUTPUT)