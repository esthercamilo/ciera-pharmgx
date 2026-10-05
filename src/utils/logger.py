import logging
import sys


def get_logger(name: str = "pipeline") -> logging.Logger:
    """Configura e retorna um logger simples com saída formatada para o console."""
    logger = logging.getLogger(name)

    # Evita duplicação de handlers caso a função seja chamada múltiplas vezes
    if not logger.handlers:
        logger.setLevel(logging.INFO)

        handler = logging.StreamHandler(sys.stdout)
        formatter = logging.Formatter(
            "%(asctime)s | %(levelname)-8s | %(name)s | %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S"
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    return logger