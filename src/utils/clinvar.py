import gzip
import xml.etree.ElementTree as ET
import csv
import time


def processar_clinvar_completo(xml_path, output_csv):
    print(f"Iniciando o processamento COMPLETO do arquivo: {xml_path}")
    print("Atenção: Por ser um arquivo grande (6GB), isso pode levar alguns minutos...")

    inicio = time.time()
    contador = 0

    # Abre o .gz em leitura e o arquivo de saída CSV em escrita simultânea
    with gzip.open(xml_path, 'rt', encoding='utf-8') as f, open(output_csv, 'w', newline='', encoding='utf-8') as out_f:
        writer = csv.writer(out_f)

        # Escreve o cabeçalho do CSV
        writer.writerow(
            ['VCV_Accession', 'Variation_Name', 'Variation_Type', 'Gene', 'Condition', 'Clinical_Significance'])

        context = ET.iterparse(f, events=('end',))

        for event, elem in context:
            if elem.tag == 'VariationArchive':

                # 1. Atributos da variante
                vcv_id = elem.get('Accession', 'N/A')
                variation_name = elem.get('VariationName', 'N/A')
                variation_type = elem.get('VariationType', 'N/A')

                # 2. Genes associados
                genes = []
                for gene in elem.findall('.//Gene'):
                    symbol = gene.get('Symbol')
                    if symbol:
                        genes.append(symbol)
                gene_str = ", ".join(set(genes)) if genes else 'N/A'

                # 3. Condições / Doenças
                conditions = []
                for element_val in elem.findall('.//Trait/Name/ElementValue'):
                    if element_val.text:
                        conditions.append(element_val.text)
                condition_str = "; ".join(set(conditions)) if conditions else 'N/A'

                # 4. Interpretação Clínica
                clnsig = elem.findtext('.//ClinicalSignificance/Description', default='N/A')

                # Escreve a linha diretamente no arquivo CSV em disco
                writer.writerow([vcv_id, variation_name, variation_type, gene_str, condition_str, clnsig])

                contador += 1

                # Mostra progresso a cada 50 mil registros para você acompanhar
                if contador % 50000 == 0:
                    print(f"Processados {contador} registros...")

                # Libera a memória do elemento XML processado
                elem.clear()

    fim = time.time()
    tempo_total = (fim - inicio) / 60
    print(f"\nProcessamento concluído com sucesso!")
    print(f"Total de registros processados: {contador:,}")
    print(f"Tempo total: {tempo_total:.2f} minutos.")
    print(f"Arquivo completo salvo em: {output_csv}")


# Caminhos dos arquivos
arquivo_xml = "/home/esther/Downloads/bioinfo/ClinVarVCVRelease_2026-10.xml.gz"
arquivo_saida = "/home/esther/Downloads/bioinfo/clinvar_vcv_completo.csv"

# Executa a extração completa
processar_clinvar_completo(arquivo_xml, arquivo_saida)