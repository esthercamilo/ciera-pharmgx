#!/usr/bin/env python3
"""
Módulo de Faseamento Farmacogenético com Beagle 5.x (Orientado a Objetos).

Painel de referência: 1000 Genomes high-coverage (CCDG, release 2020-10-28), GRCh38.
Todas as coordenadas em `DEFAULT_GENE_COORDS_GRCH38` são GRCh38/hg38. Se o VCF da
amostra estiver em GRCh37, faça liftover ANTES de chamar `phase_sample`.

Referências: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/

Principais correções em relação à versão anterior
-------------------------------------------------
* O Beagle casa marcadores por CHROM + POS + REF + ALT. A versão anterior só
  comparava (CHROM, POS); se REF/ALT da amostra diferissem do painel, o check
  dizia "1 marcador em comum" mas o Beagle via zero ("no markers in common").
  Agora há harmonização alelo-a-alelo (`_harmonize_alleles`): mantém, inverte
  (REF/ALT trocados), corrige fita (SNV) ou descarta com aviso.
* Marcadores-âncora agora têm genótipo AUSENTE (./.) em vez de 0/0. Antes o
  código "inventava" genótipos homozigotos-referência para o paciente.
* Âncoras são espalhadas pela região (antes: as 2 primeiras variantes) e só são
  usadas quando há poucos marcadores realmente compartilhados.
* Regiões que se sobrepõem (CYP2C19 x CYP2C9) geravam registros duplicados no
  concat final: coordenada do CYP2C19 ajustada e validação de sobreposição.
* Cache de download atômico (.part -> rename) e cache do recorte com a região
  no nome do arquivo; antes um download interrompido virava "cache válido".
* Entrada: só o GT é preservado (aceita FORMAT com DP:GQ etc.), multialélicos
  são separados (`bcftools norm -m -any`) e o nome real da amostra é mantido.
* `build` padrão do PharmacogeneticPhaser agora é GRCh38 (era GRCh37).
"""

from __future__ import annotations
import shutil
import subprocess
import tempfile
from pathlib import Path
from src.utils.logger import get_logger
import os

# ---------------------------------------------------------------------------
# Coordenadas GRCh38 (hg38)
# ---------------------------------------------------------------------------
# ATENÇÃO: as regiões não podem se sobrepor. O CYP2C19 terminava em 94950000,
# invadindo o início do CYP2C9 (94930000). Fim do CYP2C19 ajustado p/ 94900000
# (o gene termina em ~94.85 Mb).
DEFAULT_GENE_COORDS_GRCH38 = {
    "CYP2D6":  {"chrom": "22", "region": "22:42120000-42140000"},
    "CYP2C19": {"chrom": "10", "region": "10:94762642-94900000"},
    "CYP2C9":  {"chrom": "10", "region": "10:94930000-95000000"},
}

BASE_URL_GRCH38 = (
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
    "1000G_2504_high_coverage/working/20201028_3202_phased"
)

CHROM_FILE_TEMPLATE_GRCH38 = (
    "CCDG_14151_B01_GRM_WGS_2020-08-05_chr{chrom}."
    "filtered.shapeit2-duohmm-phased.vcf.gz"
)


# ---------------------------------------------------------------------------
# Helpers (nível de módulo, compartilhados pelas duas classes)
# ---------------------------------------------------------------------------
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def run_cmd(cmd: list[str], **kwargs) -> None:
    print(f"[cmd] {' '.join(cmd)}", flush=True)
    subprocess.run(cmd, check=True, **kwargs)


def parse_region(region: str) -> tuple[str, int, int]:
    """'chr22:100-200' -> ('chr22', 100, 200)."""
    chrom, _, rng = region.partition(":")
    start, _, end = rng.partition("-")
    return chrom, int(start.replace(",", "")), int(end.replace(",", ""))


def region_in_style(region: str, style: str) -> str:
    """Ajusta '22:100-200' -> 'chr22:100-200' (ou o inverso)."""
    chrom, sep, rest = region.partition(":")
    if style == "chr" and not chrom.startswith("chr"):
        return f"chr{chrom}{sep}{rest}"
    if style == "" and chrom.startswith("chr"):
        return f"{chrom[3:]}{sep}{rest}"
    return region


def detect_contig_style(vcf_path: Path) -> str:
    """Retorna 'chr' se o VCF usa 'chr22', senão '' (usa '22')."""
    res = subprocess.run(
        ["bcftools", "index", "-s", str(vcf_path)],
        capture_output=True, text=True,
    )
    if res.returncode == 0:
        for line in res.stdout.splitlines():
            if line.strip():
                return "chr" if line.split("\t")[0].startswith("chr") else ""
    # fallback: cabeçalho e depois primeira linha de dados
    res = subprocess.run(
        ["bcftools", "view", "-h", str(vcf_path)],
        capture_output=True, text=True, check=True,
    )
    for line in res.stdout.splitlines():
        if line.startswith("##contig") and "ID=" in line:
            cid = line.split("ID=")[1].split(",")[0].split(">")[0]
            return "chr" if cid.startswith("chr") else ""
    res = subprocess.run(
        ["bcftools", "query", "-f", "%CHROM\n", str(vcf_path)],
        capture_output=True, text=True, check=True,
    )
    first = res.stdout.split("\n", 1)[0]
    return "chr" if first.startswith("chr") else ""


def bcf_query(vcf: Path, fmt: str, region: str | None = None) -> list[list[str]]:
    """Roda `bcftools query` e devolve as linhas já separadas por TAB."""
    cmd = ["bcftools", "query"]
    if region:
        cmd += ["-r", region]
    cmd += ["-f", fmt, str(vcf)]
    out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    return [ln.split("\t") for ln in out.splitlines() if ln]


def count_variants(vcf: Path, region: str | None = None) -> int:
    return len(bcf_query(vcf, "%POS\n", region))


def _is_snv(*alleles: str) -> bool:
    return all(len(a) == 1 and a in "ACGT" for a in alleles)


def _is_palindromic(a: str, b: str) -> bool:
    """A/T e C/G: fita direta e reversa são indistinguíveis."""
    return a.translate(_COMPLEMENT) == b


def swap_gt(gt: str) -> str:
    """Inverte 0<->1 num GT bialélico ('0|1' -> '1|0', '1/1' -> '0/0')."""
    return "".join({"0": "1", "1": "0"}.get(c, c) for c in gt)


def reconcile_alleles(
    s_ref: str, s_alt: str, panel: set[tuple[str, str]]
) -> tuple[str, tuple[str, str] | None, bool]:
    """
    Compara REF/ALT da amostra com os (REF, ALT) do painel na MESMA posição.

    Retorna (status, (ref_painel, alt_painel) ou None, inverter_gt):
      exact      -> idêntico
      swap       -> REF/ALT trocados (o GT precisa ser invertido)
      flip       -> fita oposta (só SNV não palindrômico)
      flip_swap  -> fita oposta + trocado
      ambiguous  -> SNV palindrômico (A/T, C/G) trocado: não dá para distinguir
                    troca de REF/ALT de inversão de fita -> descartado
      mismatch   -> posição existe no painel, mas com outros alelos
      absent     -> posição inexistente no painel
    """
    if not panel:
        return "absent", None, False
    if (s_ref, s_alt) in panel:
        return "exact", (s_ref, s_alt), False
    if (s_alt, s_ref) in panel:
        if _is_snv(s_ref, s_alt) and _is_palindromic(s_ref, s_alt):
            return "ambiguous", None, False
        return "swap", (s_alt, s_ref), True
    if _is_snv(s_ref, s_alt):
        c_ref, c_alt = s_ref.translate(_COMPLEMENT), s_alt.translate(_COMPLEMENT)
        if (c_ref, c_alt) in panel:
            return "flip", (c_ref, c_alt), False
        if (c_alt, c_ref) in panel:
            return "flip_swap", (c_alt, c_ref), True
    return "mismatch", None, False


def write_vcf_gz(
    records: list[tuple[str, int, str, str, str, str]],
    sample_name: str,
    out_vcf: Path,
) -> None:
    """
    Escreve registros (chrom, pos, id, ref, alt, gt) em VCF.gz ordenado + índice.
    """
    contigs = sorted({r[0] for r in records})
    tmp = Path(tempfile.mktemp(suffix=".vcf", dir=str(out_vcf.parent)))
    try:
        with open(tmp, "w") as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            for c in contigs:
                f.write(f"##contig=<ID={c}>\n")
            f.write(
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                f"{sample_name}\n"
            )
            for chrom, pos, vid, ref, alt, gt in records:
                f.write(f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t.\t.\t.\tGT\t{gt}\n")
        run_cmd(["bcftools", "sort", "-Oz", "-o", str(out_vcf), str(tmp)])
        run_cmd(["bcftools", "index", "-t", "-f", str(out_vcf)])
    finally:
        tmp.unlink(missing_ok=True)


def _require_tools(*tools: str) -> None:
    for tool in tools:
        if shutil.which(tool) is None:
            raise RuntimeError(f"Ferramenta obrigatória não encontrada no PATH: {tool}")


# ---------------------------------------------------------------------------
# Painel de referência
# ---------------------------------------------------------------------------
class ReferenceManager:
    """
    Gerencia download e recorte de painéis de referência do 1000 Genomes
    para genes farmacogenéticos.

    Fluxo:
      1. Baixa (uma vez) o VCF cromossômico inteiro do 1000G GRCh38.
      2. Recorta apenas a região do gene de interesse.
      3. Indexa o VCF recortado com tabix.
      4. Reutiliza caches nas execuções seguintes.
    """

    def __init__(
        self,
        output_dir: Path | str,
        build: str = "GRCh38",
        gene_coords: dict | None = None,
    ) -> None:
        if build != "GRCh38":
            raise ValueError(
                "Esta versão do ReferenceManager só suporta GRCh38. "
                "Use um painel GRCh37 separadamente se necessário."
            )
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.build = build
        self.gene_coords = gene_coords or dict(DEFAULT_GENE_COORDS_GRCH38)

    # ------------------------------------------------------------------ utils

    def _chrom_url(self, chrom: str) -> tuple[str, str]:
        """Retorna (url_vcf, url_tbi) para um cromossomo."""
        base = chrom.replace("chr", "")
        fname = CHROM_FILE_TEMPLATE_GRCH38.format(chrom=base)
        vcf_url = f"{BASE_URL_GRCH38}/{fname}"
        return vcf_url, f"{vcf_url}.tbi"

    def _download_if_needed(self, url: str, dest: Path, label: str) -> None:
        """Baixa url -> dest se dest ainda não existir (download atômico)."""
        if dest.is_file() and dest.stat().st_size > 0:
            print(f"[info] {label} já existe em cache: {dest}")
            return
        print(f"[info] Baixando {label}...\n       → {url}\n       → {dest}")
        part = dest.with_name(dest.name + ".part")
        run_cmd(["wget", "-c", "-O", str(part), url])
        part.rename(dest)  # só vira "cache" se o download terminou

    # --------------------------------------------------------------- pipeline

    def get_gene_reference(self, gene_symbol: str) -> Path:
        """
        Retorna o path do VCF de referência recortado para o gene.
        Baixa o cromossomo inteiro na primeira execução, recorta a região
        e indexa. Nas execuções seguintes, reutiliza o cache.
        """
        _require_tools("wget", "bcftools")

        gene = gene_symbol.upper()
        if gene not in self.gene_coords:
            raise ValueError(
                f"Gene {gene} não cadastrado. "
                f"Disponíveis: {sorted(self.gene_coords.keys())}"
            )

        info = self.gene_coords[gene]
        chrom = info["chrom"]
        region = info["region"]
        _, start, end = parse_region(region)

        # A região vai no nome: se você mudar as coordenadas, o cache antigo
        # não é reaproveitado por engano.
        out_vcf = self.output_dir / (
            f"ref_{gene.lower()}_{self.build}_{start}-{end}.vcf.gz"
        )
        out_tbi = Path(str(out_vcf) + ".tbi")
        if out_vcf.is_file() and out_tbi.is_file():
            print(f"[info] Referência final para {gene} já existe em {out_vcf}")
            return out_vcf

        # 1. Baixa o cromossomo inteiro (se necessário)
        vcf_url, tbi_url = self._chrom_url(chrom)
        chrom_vcf = self.output_dir / Path(vcf_url).name
        chrom_tbi = Path(str(chrom_vcf) + ".tbi")
        self._download_if_needed(vcf_url, chrom_vcf, f"chr{chrom} VCF (GRCh38)")
        self._download_if_needed(tbi_url, chrom_tbi, f"chr{chrom} index (GRCh38)")

        # 2. Recorta a região do gene
        ref_style = detect_contig_style(chrom_vcf)
        region_style = region_in_style(region, ref_style)
        print(f"[info] Recortando região {region_style} para {gene}...")
        tmp_vcf = out_vcf.with_name("tmp_" + out_vcf.name)
        run_cmd([
            "bcftools", "view", "-r", region_style,
            "-Oz", "-o", str(tmp_vcf), str(chrom_vcf),
        ])

        # 3. Sanidade + índice
        run_cmd(["bcftools", "index", "-t", "-f", str(tmp_vcf)])
        n = count_variants(tmp_vcf)
        if n == 0:
            tmp_vcf.unlink(missing_ok=True)
            Path(str(tmp_vcf) + ".tbi").unlink(missing_ok=True)
            raise RuntimeError(
                f"O recorte para {gene} em {region} ficou vazio. "
                f"Verifique as coordenadas e o build do painel."
            )
        tmp_vcf.rename(out_vcf)
        Path(str(tmp_vcf) + ".tbi").rename(out_tbi)
        print(f"[ok] Referência GRCh38 pronta para {gene}: {out_vcf} ({n} variantes)")
        return out_vcf

    def get_all_references(self) -> dict[str, Path]:
        """Baixa e recorta todos os genes cadastrados de uma vez."""
        return {gene: self.get_gene_reference(gene) for gene in self.gene_coords}


# ---------------------------------------------------------------------------
# Faseamento
# ---------------------------------------------------------------------------
class PharmacogeneticPhaser:
    def __init__(
        self,
        beagle_jar: Path | str,
        maps_dir: Path | str,
        regions_map: dict[str, str | dict] | None = None,
        build: str = "GRCh38",
        threads: int = 4,
        memory_gb: int = 8,
        min_shared_markers: int = 2,
        n_anchors: int = 10,
    ) -> None:
        """
        min_shared_markers: se a amostra tiver menos marcadores COMPATÍVEIS com o
            painel que isso na região, são adicionados marcadores-âncora com
            genótipo ausente (imputados pelo Beagle a partir do painel).
        n_anchors: quantos marcadores-âncora adicionar nesse caso.
        """
        self.beagle_jar = Path(beagle_jar)
        self.maps_dir = Path(maps_dir)
        self.regions_map = regions_map or dict(DEFAULT_GENE_COORDS_GRCH38)
        self.build = build
        self.threads = threads
        self.memory_gb = memory_gb
        self.min_shared_markers = min_shared_markers
        self.n_anchors = n_anchors

        self._check_dependencies()

    # ------------------------------------------------------------------ utils

    def _check_dependencies(self) -> None:
        _require_tools("java", "bcftools")
        if not self.beagle_jar.is_file():
            raise FileNotFoundError(f"Beagle não encontrado em {self.beagle_jar}")

    def _find_map_file(self, chrom: str) -> Path:
        clean = chrom.replace("chr", "")
        candidates = [
            self.maps_dir / f"plink.chr{clean}.{self.build}.map",     # padrão Beagle
            self.maps_dir / f"plink.chrchr{clean}.{self.build}.map",  # coluna CHROM com 'chr'
            self.maps_dir / f"plink.{clean}.{self.build}.map",
            self.maps_dir / f"chr{clean}.{self.build}.map",
        ]
        for c in candidates:
            if c.is_file():
                return c
        raise FileNotFoundError(
            f"Mapa genético não encontrado para {chrom}. Procurei em: "
            + ", ".join(str(c) for c in candidates)
        )

    def _region_of(self, gene_or_region: str) -> str:
        val = self.regions_map[gene_or_region]
        return val["region"] if isinstance(val, dict) else val

    @staticmethod
    def _validate_no_overlap(regions: dict[str, str]) -> None:
        """Regiões sobrepostas geram registros duplicados no concat final."""
        by_chrom: dict[str, list[tuple[int, int, str]]] = {}
        for name, region in regions.items():
            chrom, start, end = parse_region(region)
            by_chrom.setdefault(chrom.replace("chr", ""), []).append((start, end, name))
        for chrom, items in by_chrom.items():
            items.sort()
            for (s1, e1, n1), (s2, e2, n2) in zip(items, items[1:]):
                if s2 <= e1:
                    raise ValueError(
                        f"Regiões sobrepostas no cromossomo {chrom}: "
                        f"{n1} ({s1}-{e1}) e {n2} ({s2}-{e2})."
                    )

    def _make_chr_map(self, add: bool, workdir: Path) -> Path:
        """Arquivo de renomeação de cromossomos para `bcftools annotate`."""
        lines = []
        for i in list(range(1, 23)) + ["X", "Y", "MT"]:
            lines.append(f"{i}\tchr{i}" if add else f"chr{i}\t{i}")
        path = workdir / ("chr_add.txt" if add else "chr_strip.txt")
        path.write_text("\n".join(lines) + "\n")
        return path

    def _harmonize_chrom(
        self, vcf_in: Path, ref_vcf: Path, workdir: Path, tag: str
    ) -> Path:
        """Reescreve o CHROM do input para casar com o estilo do painel."""
        ref_style = detect_contig_style(ref_vcf)
        in_style = detect_contig_style(vcf_in)
        if ref_style == in_style:
            return vcf_in

        print(
            f"[info] Ajustando CHROM de '{in_style or '(sem chr)'}' "
            f"para '{ref_style or '(sem chr)'}' para casar com o painel."
        )
        out = workdir / f"input_{tag}.renamed.vcf.gz"
        map_file = self._make_chr_map(add=(ref_style == "chr"), workdir=workdir)
        run_cmd([
            "bcftools", "annotate", "--rename-chrs", str(map_file),
            "-Oz", "-o", str(out), str(vcf_in),
        ])
        run_cmd(["bcftools", "index", "-t", "-f", str(out)])
        return out

    # ------------------------------------------------- harmonização de alelos

    def _harmonize_alleles(
        self, vcf_target: Path, ref_vcf: Path, region: str
    ) -> tuple[list[tuple[str, int, str, str, str, str]], int]:
        """
        Compara cada variante da amostra com o painel (CHROM, POS, REF, ALT) —
        que é o critério que o Beagle usa — e devolve os registros compatíveis
        (já com REF/ALT/GT ajustados) e quantos são.
        """
        panel: dict[tuple[str, int], set[tuple[str, str]]] = {}
        for chrom, pos, ref, alt in bcf_query(
            ref_vcf, "%CHROM\t%POS\t%REF\t%ALT\n", region
        ):
            panel.setdefault((chrom, int(pos)), set()).add((ref, alt))

        kept: list[tuple[str, int, str, str, str, str]] = []
        for chrom, pos, vid, ref, alt, gt in bcf_query(
            vcf_target, "%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n"
        ):
            pos = int(pos)
            here = panel.get((chrom, pos), set())
            status, pair, invert = reconcile_alleles(ref, alt, here)
            label = f"{chrom}:{pos} {ref}>{alt} (GT {gt})"
            if pair is None:
                painel_txt = ", ".join(f"{r}>{a}" for r, a in sorted(here)) or "—"
                print(
                    f"[descartado] {label}: {status}. "
                    f"Alelos do painel nessa posição: {painel_txt}"
                )
                continue
            new_gt = swap_gt(gt) if invert else gt
            if status != "exact":
                print(
                    f"[ajustado] {label}: {status} -> "
                    f"{pair[0]}>{pair[1]} (GT {new_gt})"
                )
            kept.append((chrom, pos, vid, pair[0], pair[1], new_gt))
        return kept, len(kept)

    def _anchor_records(
        self,
        ref_vcf: Path,
        region: str,
        sample_positions: set[int],
        n: int,
    ) -> list[tuple[str, int, str, str, str, str]]:
        """
        Escolhe `n` SNVs do painel, espalhados pela região, e cria registros com
        genótipo AUSENTE (./.). Genótipo ausente não carrega informação: serve só
        para o Beagle ter marcadores em comum e imputar a partir do painel.
        """
        cands = [
            (chrom, int(pos), ref, alt)
            for chrom, pos, ref, alt in bcf_query(
                ref_vcf, "%CHROM\t%POS\t%REF\t%ALT\n", region
            )
            if _is_snv(ref, alt) and int(pos) not in sample_positions
        ]
        if not cands:
            raise RuntimeError(
                f"Nenhum SNV do painel disponível em {region} para servir de âncora."
            )
        n = min(n, len(cands))
        idx = sorted({round(i * (len(cands) - 1) / max(n - 1, 1)) for i in range(n)})
        return [
            (cands[i][0], cands[i][1], ".", cands[i][2], cands[i][3], "./.")
            for i in idx
        ]

    # --------------------------------------------------------------- pipeline

    def prepare_input(self, vcf_in: Path | str) -> Path:
        """
        Lê o VCF da amostra (tolera separador por espaços), mantém só
        CHROM/POS/ID/REF/ALT + GT, ordena, separa multialélicos e indexa.
        """
        vcf_path = Path(vcf_in)
        workdir = vcf_path.parent / "output_work"
        workdir.mkdir(parents=True, exist_ok=True)

        sample_name = "Sample"
        records: list[tuple[str, int, str, str, str, str]] = []

        with open(vcf_path, "r") as fin:
            for raw in fin:
                line = raw.rstrip("\r\n")
                if not line.strip() or line.startswith("##"):
                    continue
                if line.startswith("#"):
                    cols = line.split()
                    if len(cols) > 9:
                        sample_name = cols[9]
                    if len(cols) > 10:
                        print(
                            f"[aviso] VCF com {len(cols) - 9} amostras; "
                            f"usando apenas a primeira ({sample_name})."
                        )
                    continue
                cols = line.split("\t") if "\t" in line else line.split()
                if len(cols) < 5:
                    raise ValueError(f"Linha de VCF inválida (<5 colunas): {line!r}")
                chrom, pos, vid, ref, alt = cols[:5]
                gt = "./."
                if len(cols) >= 10:
                    fmt = cols[8].split(":")
                    if "GT" in fmt:
                        gt = cols[9].split(":")[fmt.index("GT")]
                else:
                    print(f"[aviso] {chrom}:{pos} sem coluna de genótipo; GT ausente.")
                records.append((chrom, int(pos), vid, ref, alt, gt))

        if not records:
            raise ValueError(f"Nenhuma variante encontrada em {vcf_path}")

        clean_vcf = workdir / f"{vcf_path.stem}.clean.vcf.gz"
        print(f"[info] Ordenando e compactando: {clean_vcf}")
        write_vcf_gz(records, sample_name, clean_vcf)

        # Separa multialélicos (Beagle/painel são bialélicos)
        normed = workdir / f"{vcf_path.stem}.sorted.vcf.gz"
        run_cmd(["bcftools", "norm", "-m", "-any", "-Oz", "-o", str(normed), str(clean_vcf)])
        run_cmd(["bcftools", "index", "-t", "-f", str(normed)])
        return normed

    def phase_sample(
        self,
        vcf_in: Path | str,
        out_prefix: Path | str,
        ref_panels: dict[str, Path | str],
        keep_workdir: bool = False,
    ) -> Path:
        if not isinstance(ref_panels, dict) or not ref_panels:
            raise ValueError("'ref_panels' deve ser um dicionário gene -> VCF.")
        unknown = [g for g in ref_panels if g not in self.regions_map]
        if unknown:
            raise ValueError(f"Região genômica não mapeada para o(s) alvo(s): {unknown}")
        self._validate_no_overlap({g: self._region_of(g) for g in ref_panels})

        vcf_in = self.prepare_input(vcf_in)
        out_prefix = Path(out_prefix)
        out_prefix.parent.mkdir(parents=True, exist_ok=True)
        workdir = Path(f"{out_prefix}_work")
        workdir.mkdir(parents=True, exist_ok=True)

        sample_name = subprocess.run(
            ["bcftools", "query", "-l", str(vcf_in)],
            capture_output=True, text=True, check=True,
        ).stdout.split()[0]

        phased_vcfs: list[Path] = []

        for gene_or_region, ref_vcf in ref_panels.items():
            print(f"\n=== Processando Alvo: {gene_or_region} ===")

            ref_vcf_path = Path(ref_vcf)
            if not ref_vcf_path.is_file():
                raise FileNotFoundError(
                    f"Arquivo de referência para {gene_or_region} não encontrado: {ref_vcf_path}"
                )
            region_grch = self._region_of(gene_or_region)

            # 1) Harmoniza CHROM do input com o do painel e ajusta a região
            vcf_in_h = self._harmonize_chrom(vcf_in, ref_vcf_path, workdir, gene_or_region)
            ref_style = detect_contig_style(ref_vcf_path)
            region_style = region_in_style(region_grch, ref_style)

            # 2) Extrai a região da amostra
            print(f"[info] Extraindo região {region_style} do VCF da amostra...")
            vcf_region = workdir / f"input_{gene_or_region}.region.vcf.gz"
            run_cmd([
                "bcftools", "view", "-r", region_style,
                "-Oz", "-o", str(vcf_region), str(vcf_in_h),
            ])
            run_cmd(["bcftools", "index", "-t", "-f", str(vcf_region)])

            n_input = count_variants(vcf_region)
            n_ref = count_variants(ref_vcf_path, region_style)
            print(f"[info] Variantes da amostra na região: {n_input}")
            print(f"[info] Variantes do painel na região:  {n_ref}")
            if n_ref == 0:
                raise RuntimeError(
                    f"O painel {ref_vcf_path} não tem variantes em {region_style}. "
                    f"Confira build (GRCh37 vs GRCh38) e coordenadas."
                )

            # 3) Casa (CHROM, POS, REF, ALT) com o painel — o critério do Beagle
            records, n_shared = self._harmonize_alleles(
                vcf_region, ref_vcf_path, region_style
            )
            print(f"[info] Marcadores compatíveis com o painel (pos+REF+ALT): {n_shared}")

            if n_shared == 0:
                print(
                    "[AVISO] NENHUM genótipo observado da amostra é compatível com o "
                    "painel nesta região. O resultado será apenas imputação a partir "
                    "do painel (sem dado real do paciente): não use para chamar "
                    "haplótipos/alelos-estrela."
                )

            # 4) Poucos marcadores compatíveis -> âncoras com genótipo ausente
            if n_shared < self.min_shared_markers:
                print(
                    f"[aviso] Só {n_shared} marcador(es) compatível(is) (mínimo "
                    f"{self.min_shared_markers}). Adicionando {self.n_anchors} "
                    f"âncoras do painel com genótipo ausente (./.); o Beagle "
                    f"imputa esses sítios a partir do painel."
                )
                anchors = self._anchor_records(
                    ref_vcf_path, region_style,
                    sample_positions={r[1] for r in records},
                    n=self.n_anchors,
                )
                records = records + anchors

            vcf_target = workdir / f"input_{gene_or_region}.vcf.gz"
            write_vcf_gz(records, sample_name, vcf_target)

            # 5) Beagle
            base_chrom = region_grch.split(":")[0].replace("chr", "")
            map_file = self._find_map_file(base_chrom)
            out_target_prefix = workdir / f"phased_{gene_or_region}"
            cmd = [
                "java", f"-Xmx{self.memory_gb}g", "-jar", str(self.beagle_jar),
                f"gt={vcf_target}",
                f"ref={ref_vcf_path}",
                f"map={map_file}",
                f"chrom={region_style}",
                f"out={out_target_prefix}",
                f"nthreads={self.threads}",
            ]
            print(f"[info] Painel de referência: {ref_vcf_path}")
            run_cmd(cmd)

            phased_vcf = Path(str(out_target_prefix) + ".vcf.gz")
            if not phased_vcf.is_file():
                raise RuntimeError(f"Beagle falhou ao gerar {phased_vcf} para {gene_or_region}")
            run_cmd(["bcftools", "index", "-t", "-f", str(phased_vcf)])
            phased_vcfs.append(phased_vcf)

        out_final = Path(str(out_prefix) + ".pharmacogenes_phased.vcf.gz")
        if len(phased_vcfs) == 1:
            shutil.copy(phased_vcfs[0], out_final)
        else:
            run_cmd([
                "bcftools", "concat", "-a", "-Oz", "-o", str(out_final),
                *[str(v) for v in phased_vcfs],
            ])
        run_cmd(["bcftools", "index", "-t", "-f", str(out_final)])
        print(f"\n[ok] VCF farmacogenético faseado gerado com sucesso: {out_final}")

        if not keep_workdir:
            shutil.rmtree(workdir, ignore_errors=True)
            print(f"[info] Diretório temporário limpo: {workdir}")

        return out_final

logger = get_logger(__name__)


def run_phasing_pipeline(vcf_in: str | Path, out_prefix: str | Path) -> Path:
    """
    Executa o faseamento genômico. Se a variável de ambiente MOCK=1 estiver ativa,
    realiza apenas a troca simulada de '/' para '|' (modo leve para testes rápidos).
    Caso contrário, executa o pipeline completo com Beagle e 1000 Genomes.
    """
    vcf_path = Path(vcf_in)
    out_prefix_path = Path(out_prefix)
    out_prefix_path.parent.mkdir(parents=True, exist_ok=True)

    # Verifica se o modo mock/simples via variável de ambiente está ativado
    is_mock = os.getenv("MOCK", "0") == "1"

    if is_mock:
        logger.warning(
            "[MOCK ATIVO] Variável MOCK=1 detectada. Pulando o Beagle e realizando apenas a normalização de barras ('/' -> '|')...")

        out_final = Path(str(out_prefix_path) + ".pharmacogenes_phased.vcf.gz")

        # Leitura simples do VCF de entrada e reescrita trocando o separador de genótipo
        import gzip
        import shutil

        temp_out = out_prefix_path.with_suffix(".mock_phased.vcf")
        with open(vcf_path, "r", encoding="utf-8", errors="ignore") as fin, open(temp_out, "w",
                                                                                 encoding="utf-8") as fout:
            for line in fin:
                if line.startswith("#"):
                    fout.write(line)
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 10:
                    # Ajusta o formato do genótipo trocando '/' por '|' se presente
                    fmt_cols = parts[8].split(":")
                    if "GT" in fmt_cols:
                        gt_idx = fmt_cols.index("GT")
                        sample_vals = parts[9].split(":")
                        gt_val = sample_vals[gt_idx]
                        if "/" in gt_val:
                            sample_vals[gt_idx] = gt_val.replace("/", "|")
                            parts[9] = ":".join(sample_vals)
                fout.write("\t".join(parts) + "\n")

        # Compacta e indexa com bcftools
        import subprocess
        subprocess.run(["bcftools", "sort", "-Oz", "-o", str(out_final), str(temp_out)], check=True)
        subprocess.run(["bcftools", "index", "-t", "-f", str(out_final)], check=True)
        temp_out.unlink(missing_ok=True)

        logger.info(f"[MOCK ATIVO] VCF pseudo-faseado gerado com sucesso: {out_final}")
        return out_final

    # --- Fluxo de Produção Real (Beagle + 1000G) ---
    logger.info("Iniciando pipeline de faseamento real (Beagle 5.x)...")

    ref_output_dir = "/media/esther/SAMSUNG/1000_genomes"
    beagle_jar_path = "/home/esther/GitProjects/flexlab/backend/genomics/softwares/beagle.27Feb25.75f.jar"
    maps_dir_path = "/home/esther/GitProjects/flexlab/backend/genomics/data/chr_in_chrom_field"

    manager = ReferenceManager(
        output_dir=ref_output_dir,
        build="GRCh38",
    )

    references = {
        "CYP2D6": manager.get_gene_reference("CYP2D6"),
        "CYP2C19": manager.get_gene_reference("CYP2C19"),
        "CYP2C9": manager.get_gene_reference("CYP2C9"),
    }

    phaser = PharmacogeneticPhaser(
        beagle_jar=beagle_jar_path,
        maps_dir=maps_dir_path,
        regions_map=manager.gene_coords,
        build="GRCh38",
    )

    phased_vcf = phaser.phase_sample(
        vcf_in=vcf_path,
        out_prefix=out_prefix_path,
        ref_panels=references,
    )

    logger.info(f"Faseamento real concluído com sucesso. VCF gerado: {phased_vcf}")
    return phased_vcf

if __name__ == "__main__":
    manager = ReferenceManager(
        "/media/esther/SAMSUNG/1000_genomes",
        build="GRCh38",
    )
    references = {
        "CYP2D6": manager.get_gene_reference("CYP2D6"),
        "CYP2C19": manager.get_gene_reference("CYP2C19"),
        "CYP2C9": manager.get_gene_reference("CYP2C9"),
    }

    phaser = PharmacogeneticPhaser(
        beagle_jar=(
            "/home/esther/GitProjects/flexlab/backend/genomics/softwares/"
            "beagle.27Feb25.75f.jar"
        ),
        maps_dir=(
            "/home/esther/GitProjects/flexlab/backend/genomics/data/"
            "chr_in_chrom_field"
        ),
        regions_map=manager.gene_coords,
        build="GRCh38",
    )

    phaser.phase_sample(
        vcf_in=(
            "/home/esther/GitProjects/flexlab/backend/genomics/tests/"
            "f32_test1/file2.vcf"
        ),
        out_prefix=(
            "/home/esther/GitProjects/flexlab/backend/genomics/tests/"
            "f32_test1/output"
        ),
        ref_panels=references,
    )