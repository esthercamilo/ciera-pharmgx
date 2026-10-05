import argparse
import json
import sys
from utils.logger import get_logger

logger = get_logger(__name__)


def parse_arguments():
    """Configura e lê os argumentos passados via linha de comando."""
    parser = argparse.ArgumentParser(
        description="Pipeline Farmacogenético - Análise e Imputação"
    )
    parser.add_argument(
        "--vcf",
        required=True,
        help="Caminho para o arquivo VCF de entrada (ex: data/input/amostra.vcf)"
    )
    parser.add_argument(
        "--patient",
        required=True,
        help="Caminho para o arquivo JSON com os metadados da amostra do paciente (ex: data/input/patient.json)"
    )
    parser.add_argument(
        "-o", "--output",
        default="data/output/resultado.json",
        help="Caminho para o arquivo JSON de saída (padrão: data/output/resultado.json)"
    )
    return parser.parse_args()


def run_mock_pipeline(vcf_path: str, patient_path: str) -> dict:
    """Simula a execução do pipeline integrando VCF e dados clínicos do paciente."""
    logger.info(f"Processando arquivo VCF: {vcf_path}")
    logger.info(f"Lendo dados do paciente: {patient_path}")

    # Exemplo de leitura opcional dos dados do paciente para enriquecer o mock
    try:
        with open(patient_path, "r", encoding="utf-8") as f:
            patient_data = json.load(f)
    except Exception as e:
        logger.warning(f"Não foi possível ler o JSON do paciente ({e}), prosseguindo com dados padrão.")
        patient_data = {"patient_id": "unknown"}

    # Mock do resultado farmacogenético consolidado
    mock_result = {
        "status": "success",
        "metadata": {
            "vcf_file": vcf_path,
            "patient_info": patient_data
        },
        "variants_analyzed": 42,
        "star_alleles": {
            "CYP2D6": "*1/*4",
            "CYP2C19": "*1/*2",
            "VKORC1": "-1639G>A"
        },
        "phenotypes": {
            "CYP2D6": "Intermediate Metabolizer",
            "CYP2C19": "Intermediate Metabolizer",
            "VKORC1": "Normal Sensitivity"
        },
        "recommendations": [
            "Considere ajuste de dose para substratos de CYP2D6 conforme diretrizes CPIC."
        ]
    }

    return mock_result


def main():
    args = parse_arguments()
    logger.info("Pipeline iniciado.")

    try:
        # Executa a lógica integrando VCF e dados do paciente
        resultado = run_mock_pipeline(args.vcf, args.patient)

        # Salva o output em JSON
        with open(args.output, "w", encoding="utf-8") as f:
            json.dump(resultado, f, indent=4, ensure_ascii=False)

        logger.info(f"Pipeline finalizado com sucesso. Output salvo em: {args.output}")

    except Exception as e:
        logger.error(f"Erro crítico durante a execução do pipeline: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()