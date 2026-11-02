# Pipeline Farmacogenético Containerizado

Pipeline para análise de variantes farmacogenéticas, integrando arquivos VCF e metadados clínicos do paciente, utilizando Python, Docker e Beagle.

---

## Pré-requisitos

Certifique-se de ter as seguintes ferramentas instaladas na sua máquina host:
* [Docker](https://docs.docker.com/get-docker/)
* [Docker Compose](https://docker.com/) (ou o plugin moderno `docker compose`)

---

## Estrutura do Projeto

```text
.
├── docker-compose.yml
├── Dockerfile
├── README.md
├── external_tools/
│   └── beagle/
│       └── beagle.jar          # Coloque o executável do Beagle aqui
├── data/                       
│   ├── input/                  # Coloque aqui seus .vcf e o .json do paciente
│   └── output/                 # Onde os resultados JSON serão salvos
└── src/                        
    ├── main.py                 # Orquestrador principal do pipeline
    ├── config.py               # Configurações globais
    ├── modules/                # Módulos lógicos (preprocessor, annotator, etc.)
    └── utils/                  
        └── logger.py           # Logger simples
        
```

## Configuração Inicial

1. **Configuração de Variáveis de Ambiente:**
   Crie um arquivo `.env` na raiz do projeto com as configurações necessárias para a execução do pipeline (por exemplo, caminhos ou parâmetros globais).

2. **Dados de Entrada:**
   O projeto já inclui dados de teste na pasta `data/input/`. Caso deseje analisar novas amostras, basta colocar os novos arquivos VCF e o JSON de metadados do paciente no mesmo diretório `data/input/`.


## Como Fazer o Build

Para construir a imagem Docker do pipeline pela primeira vez (ou sempre que atualizar dependências):

```bash
$ sudo docker compose build
```

*(Nota: Se o seu usuário estiver no grupo `docker`, você pode omitir o `sudo`).*


## Como Executar (Amostra Individual)

Para rodar o pipeline para uma amostra específica via linha de comando, utilizando um container temporário (`--rm`):

```bash
$ sudo docker compose run --rm pharmacogenomics-pipeline python src/main.py \
    --vcf data/input/amostra.vcf \
    --patient data/input/patient.json \
    -o data/output/resultado.json
```

O resultado processado estará disponível no host em `data/output/resultado.json`.

## Nota sobre Escala e Processamento em Lote

O comando `docker compose run` é excelente para validações pontuais, testes de desenvolvimento e processamento de amostras avulsas. No entanto, **ele não é viável para cenários de produção com centenas ou milhares de pacientes**, pois exigiría a execução manual ou scripts externos desnecessariamente complexos para cada arquivo.

### Como essa limitação será mitigada:
Para suportar alto volume de dados, a arquitetura evoluirá para um **script orquestrador em lote (Batch Runner)** executado dentro do container. Esse script será responsável por:
1. Varrer a pasta `data/input/` automaticamente.
2. Identificar e parear todos os arquivos VCFs e JSONs de pacientes.
3. Executar o pipeline de forma automatizada (utilizando paralelismo ou filas) iterando por toda a base de dados em uma única chamada.

## Comando para restore do banco

```
psql "postgresql://pregenescan:password@dpg-db9m7eom7kps73dtm8t0-a.ohio-postgres.render.com/ciera_pgx_db?sslmode=require" -c "DROP SCHEMA public CASCADE; CREATE SCHEMA public;" && psql "postgresql://pregenescan:2Tqm9XnOtk8QMWI2oljq7fnp3jggIvLZ@dpg-db9m7eom7kps73dtm8t0-a.ohio-postgres.render.com/ciera_pgx_db?sslmode=require" -f ciera_pgx_db_render_2026-11-01.dump
```

## Comando para backup do banco

```
pg_dump -h dpg-db9m7eom7kps73dtm8t0-a.ohio-postgres.render.com -U pregenescan ciera_pgx_db > ciera_pgx_db_render_2026-11-01b.dump
```