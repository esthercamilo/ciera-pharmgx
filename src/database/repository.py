import os
import psycopg2
from dotenv import load_dotenv
from src.utils.logger import get_logger

# Carrega as variáveis do arquivo .env se ele existir na raiz
load_dotenv()

logger = get_logger(__name__)


class DatabaseConnection:
    """Gerenciador de contexto para gerenciar a conexão com o PostgreSQL usando variáveis do .env."""

    def __init__(self):
        self.conn = None

    def __enter__(self):
        self.conn = psycopg2.connect(
            dbname=os.getenv("DB_NAME"),
            user=os.getenv("DB_USER"),
            password=os.getenv("DB_PASSWORD"),
            host=os.getenv("DB_HOST"),
            port=os.getenv("DB_PORT", "5432"),
            sslmode=os.getenv("DB_SSLMODE", "require"),
        )
        return self.conn

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.conn:
            self.conn.close()


class PharmacogeneticRepository:
    """Repositório especializado para consultas farmacogenéticas."""

    @staticmethod
    def get_diplotype_phenotype(conn, gene: str, diplotype_str: str) -> dict:
        """Busca o score e fenótipo de um diplótipo para um dado gene."""
        cursor = conn.cursor()
        try:
            query = """
                SELECT d.score, d.descricao 
                FROM diplotipos d
                JOIN haplotipos h ON d.hap1 = h.id
                JOIN genes g ON h.gene_id = g.id
                WHERE UPPER(g.nome) = UPPER(%s) AND d.descricao ILIKE %s
                LIMIT 1;
            """
            cursor.execute(query, (gene, f"%{diplotype_str}%"))
            row = cursor.fetchone()

            if row:
                return {
                    "score": float(row[0]) if row[0] is not None else None,
                    "phenotype": row[1]
                }
            return {"score": None, "phenotype": "Desconhecido / Não catalogado"}
        finally:
            cursor.close()