import os
import shutil
from jinja2 import Template
from src.utils.logger import get_logger

logger = get_logger(__name__)


def generate_html_preview(output_json_path: str, report_data: dict) -> str:
    """Renderiza o HTML usando Jinja2 e copia o CSS e o Logo para o diretório de output."""
    html_path = output_json_path.replace(".json", ".html")
    output_dir = os.path.dirname(output_json_path)
    logger.info(f"Gerando arquivo HTML em: {html_path}")

    # Descobre a raiz do projeto (sobe duas pastas de src/utils)
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    reports_dir = os.path.join(base_dir, "src", "reports")

    template_path = os.path.join(reports_dir, "template1.html")
    css_source_path = os.path.join(reports_dir, "style.css")
    logo_source_path = os.path.join(reports_dir, "logo_ciera.avif")

    if not os.path.exists(template_path):
        raise FileNotFoundError(f"Template HTML não encontrado em: {template_path}")

    # Copia o style.css para o diretório de output
    if os.path.exists(css_source_path):
        shutil.copyfile(css_source_path, os.path.join(output_dir, "style.css"))
        logger.info("Arquivo style.css copiado com sucesso.")

    # Copia o logo_ciera.avif para o diretório de output
    if os.path.exists(logo_source_path):
        shutil.copyfile(logo_source_path, os.path.join(output_dir, "logo_ciera.avif"))
        logger.info("Arquivo logo_ciera.avif copiado com sucesso.")
    else:
        logger.warning(f"Logo não encontrado em {logo_source_path}. A imagem pode ficar quebrada.")

    # Lê o template HTML
    with open(template_path, "r", encoding="utf-8") as f:
        template_content = f.read()

    # Compila o Jinja2
    template = Template(template_content)

    # Renderiza passando o dicionário completo de contexto
    rendered_html = template.render(
        patient=report_data.get("patient_summary", {}),
        requester=report_data.get("requester_info", {}),
        dates=report_data.get("dates", {}),
        summary_impact=report_data.get("summary_impact", []),
        detailed_findings=report_data.get("detailed_findings", []),
        references=report_data.get("references", []),
        signature=report_data.get("signature", {})
    )

    # Escreve o arquivo HTML final
    with open(html_path, "w", encoding="utf-8") as f:
        f.write(rendered_html)

    logger.info(f"HTML gerado com sucesso: {html_path}")
    return html_path