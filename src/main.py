import argparse
import json
import os
import sys
from pathlib import Path

from src.modules.parse_vcf import enrich_variants_from_database, parse_vcf_variants
from src.reports.html_renderer import generate_html_preview
from src.utils.logger import get_logger

# Importação da função de faseamento que criamos (com suporte a MOCK=1)
# from src.modules.phasing_service import run_phasing_pipeline

logger = get_logger(__name__)


def parse_arguments():
    """Configura e lê os argumentos passados via linha de comando."""
    parser = argparse.ArgumentParser(
        description="Pipeline Ciera-PharmGX - Análise Farmacogenômica Real"
    )
    parser.add_argument(
        "--vcf",
        required=True,
        help="Caminho para o arquivo VCF de entrada"
    )
    parser.add_argument(
        "--patient",
        required=True,
        help="Caminho para o arquivo JSON com os metadados do paciente"
    )
    parser.add_argument(
        "-o", "--output-dir",
        default=None,
        help="Diretório de saída para os resultados. Se omitido, usa a mesma pasta do VCF de entrada."
    )
    parser.add_argument(
        "--preview-html",
        action="store_true",
        help="Indica se deve gerar a pré-visualização em HTML junto com o JSON"
    )
    return parser.parse_args()


def build_report_payload(patient_data: dict, clinical_enrichment: dict, references: list) -> dict:
    """
    Constrói o payload do relatório unindo os dados do paciente
    e o enriquecimento clínico gerado pelo pipeline, sem restrições ou dados fixos.
    """
    return {
        "patient": patient_data,
        "clinical_findings": clinical_enrichment.get("detailed_findings", []),
        "summary_impact": clinical_enrichment.get("summary_impact", []),
        "references": references
    }


def main():
    args = parse_arguments()
    logger.info("Iniciando pipeline Ciera-PharmGX...")

    try:
        base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

        # 1. Resolução dos caminhos de entrada
        vcf_path = args.vcf if os.path.isabs(args.vcf) else os.path.join(base_dir, args.vcf)
        patient_path = args.patient if os.path.isabs(args.patient) else os.path.join(base_dir, args.patient)

        if not os.path.exists(vcf_path):
            raise FileNotFoundError(f"Arquivo VCF não encontrado em: {vcf_path}")
        if not os.path.exists(patient_path):
            raise FileNotFoundError(f"Arquivo de metadados do paciente não encontrado em: {patient_path}")

        # 2. Definição inteligente do diretório de saída
        if args.output_dir:
            output_dir = args.output_dir if os.path.isabs(args.output_dir) else os.path.join(base_dir, args.output_dir)
        else:
            # Se não informado, usa a mesma pasta do VCF de entrada criando uma subpasta 'output' (ou na raiz do input)
            input_folder = os.path.dirname(vcf_path)
            output_dir = os.path.join(input_folder, "output")

        os.makedirs(output_dir, exist_ok=True)
        logger.info(f"Diretório de saída definido em: {output_dir}")

        # Define os caminhos dos arquivos finais com base no nome do VCF de entrada
        vcf_stem = Path(vcf_path).stem.replace(".vcf", "")
        output_json_path = os.path.join(output_dir, f"{vcf_stem}_resultado.json")

        # 3. Leitura dos Metadados do Paciente
        with open(patient_path, "r", encoding="utf-8") as f:
            patient_data = json.load(f)
        logger.info(f"Metadados do paciente carregados de: {patient_path}")

        # 4. Etapa de Faseamento (Suporta MOCK=1 ou execução real)
        logger.info("Iniciando etapa de faseamento genômico...")
        # phased_vcf_path = run_phasing_pipeline(vcf_path, os.path.join(output_dir, vcf_stem))
        # Para testes atuais, caso prefira apontar direto para o VCF se já estiver faseado:
        phased_vcf_path = vcf_path

        # 5. Parsing do VCF Faseado
        logger.info(f"Extraindo variantes do VCF: {phased_vcf_path}")
        raw_variants = parse_vcf_variants(phased_vcf_path)

        # 6. Consulta ao Banco de Dados (Enriquecimento Clínico)
        logger.info("Realizando cruzamento clínico com o banco de dados...")
        clinical_enrichment = enrich_variants_from_database(raw_variants)
        references = []  # todo

        # 7. Montagem do Laudo
        report_payload = build_report_payload(patient_data, clinical_enrichment, references)

        # 8. Persistência do JSON de Saída
        with open(output_json_path, "w", encoding="utf-8") as f:
            json.dump(report_payload, f, indent=4, ensure_ascii=False)
        logger.info(f"JSON gerado com sucesso em: {output_json_path}")

        # 9. Renderização Opcional do HTML
        if args.preview_html:
            logger.info("Gerando pré-visualização em HTML...")
            generate_html_preview(output_json_path, report_payload)

        logger.info("Pipeline executado com sucesso.")

    except Exception as e:
        logger.error(f"Erro crítico na execução do pipeline: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()