FROM python:3.10-slim

# Instalar dependências do sistema, incluindo Java para o Beagle
RUN apt-get update && apt-get install -y --no-install-recommends \
    default-jre \
    procps \
    && rm -rf /var/lib/apt/lists/*

# Definir diretório de trabalho
WORKDIR /app

# Instalar dependências Python do projeto (se houver)
COPY src/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copiar código fonte
COPY src/ /app/src/

# Criar diretórios de dados padrão
RUN mkdir -p /app/data/input /app/data/output /app/external_tools/beagle

# Definir variáveis de ambiente úteis
ENV PYTHONUNBUFFERED=1

CMD ["python", "src/pipeline.py"]