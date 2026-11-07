#!/bin/bash

# Verifica se os dois argumentos foram passados
if [ "$#" -ne 2 ]; then
    echo "Uso: ./ciera-run.sh <caminho_para_vcf> <caminho_para_paciente_json>"
    exit 1
fi

# Pega os caminhos absolutos dos arquivos de entrada
VCF_PATH=$(realpath "$1")
PATIENT_PATH=$(realpath "$2")

# Diretórios base de entrada
VCF_DIR=$(dirname "$VCF_PATH")
PATIENT_DIR=$(dirname "$PATIENT_PATH")

# Sua pasta externa definitiva para dados (fora do projeto)
EXTERNAL_DATA_DIR="/home/esther/Documents/data"

# Cria a pasta de saída externa caso ela não exista
mkdir -p "$EXTERNAL_DATA_DIR/output"

# Executa o container mapeando:
# 1. As pastas dos arquivos de entrada
# 2. A pasta externa inteira em /app/data para que input e output fiquem sincronizados
sudo docker run --rm \
  --env-file .env \
  -v "$VCF_DIR:/data/input_vcf" \
  -v "$PATIENT_DIR:/data/input_patient" \
  -v "$EXTERNAL_DATA_DIR:/app/data" \
  ciera-pharmgx:latest \
  env PYTHONPATH=/app python src/main.py \
    --vcf "/data/input_vcf/$(basename "$VCF_PATH")" \
    --patient "/data/input_patient/$(basename "$PATIENT_PATH")" \
    --output /app/data/output/resultado.json \
    --preview-html

sudo chown -R $(id -u):$(id -g) "$EXTERNAL_DATA_DIR/output"

echo "Processamento concluído! Verifique em: $EXTERNAL_DATA_DIR/output/"