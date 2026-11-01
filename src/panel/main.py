import os
import re

def parse_bed_interval(line):
    parts = line.strip().split('\t')
    if len(parts) < 3:
        return None
    chrom = parts[0]
    try:
        start = int(parts[1])
        end = int(parts[2])
    except ValueError:
        return None
    name = parts[3] if len(parts) > 3 else "."
    extra = "\t".join(parts[4:]) if len(parts) > 4 else "."
    return chrom, start, end, name, extra


def overlaps(c1, s1, e1, c2, s2, e2):
    if c1 != c2:
        return False
    return not (e1 <= s2 or s1 >= e2)


def clean_to_single_id(name):
    """Garante que o rsID/nome seja curto e válido para o AmpliSeq (sem vírgulas ou caracteres inválidos)."""
    if not name or name == ".":
        return "target"
    # Pega apenas o primeiro ID antes de vírgula ou barra
    first_item = re.split(r'[,;/]+', str(name))[0].strip()
    # Remove tudo que não for alfanumérico, hífen ou underline
    sanitized = re.sub(r'[^a-zA-Z0-9_\-]', '', first_item)
    return sanitized if sanitized else "target"


def generate_cr_pgx_v9():
    missed_file = "IAA41053_181_Missed.bed"
    designed_file = "IAA41053_181_Designed.bed"
    illumina_probes_file = "3647_hg38_PGx_Probes_229705.bed"
    output_file = "cr-pgx-v9.bed"

    if not os.path.exists(illumina_probes_file):
        print(f"[ERRO CRÍTICO] O arquivo de probes da Illumina '{illumina_probes_file}' não foi encontrado na pasta!")
        return

    # 1. Carregar as regiões que falharam (Missed)
    missed_regions = []
    if os.path.exists(missed_file):
        with open(missed_file, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parsed = parse_bed_interval(line)
                if parsed:
                    missed_regions.append(parsed)

    # 2. Ler Designed e remover os que estão nos Missed (por posição)
    valid_entries = []
    removed_entries = []
    id_counter = 589352

    if os.path.exists(designed_file):
        with open(designed_file, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    extra_field = parts[5] if len(parts) > 5 else ""

                    rs_id = "."
                    if "SUBMITTED_REGION=" in extra_field:
                        for chunk in extra_field.split(';'):
                            if chunk.startswith("SUBMITTED_REGION="):
                                rs_id = chunk.split('=')[1]
                    elif len(parts) > 4:
                        rs_id = parts[4]

                    is_missed = False
                    for m_chr, m_start, m_end, m_name, m_extra in missed_regions:
                        if overlaps(chrom, start, end, m_chr, m_start, m_end):
                            is_missed = True
                            break

                    cleaned_rs_id = clean_to_single_id(rs_id)

                    if is_missed:
                        removed_entries.append((chrom, start, end, cleaned_rs_id))
                    else:
                        valid_entries.append((chrom, start, end, str(id_counter), cleaned_rs_id, "."))
                        id_counter += 1

    # 3. Ler EXCLUSIVAMENTE o arquivo de Probes da Illumina para pegar os amplicons de ~120pb corretos
    added_entries = []
    with open(illumina_probes_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parsed = parse_bed_interval(line)
            if parsed:
                il_chrom, il_start, il_end, il_name, il_extra = parsed

                # Verificar se o probe da Illumina intercepta alguma região do Missed
                should_add = False
                matched_rs = "."
                for m_chr, m_start, m_end, m_name, m_extra in missed_regions:
                    if overlaps(il_chrom, il_start, il_end, m_chr, m_start, m_end):
                        should_add = True
                        matched_rs = m_extra if "_" in m_extra else (m_name if "_" in m_name else il_name)
                        break

                if should_add:
                    # Garantir que o probe tem tamanho realista (ex: entre 50 e 300 pb)
                    size = il_end - il_start
                    if size <= 300:
                        cleaned_matched_rs = clean_to_single_id(matched_rs)
                        new_entry = (il_chrom, il_start, il_end, str(id_counter), cleaned_matched_rs, ".")
                        if new_entry not in added_entries:
                            added_entries.append(new_entry)
                            id_counter += 1

    final_entries = valid_entries + added_entries

    def chrom_sort_key(item):
        c = item[0]
        c_num = c.replace('chr', '')
        if c_num.isdigit():
            return (0, int(c_num), item[1])
        elif c_num == 'X':
            return (1, 23, item[1])
        elif c_num == 'Y':
            return (2, 24, item[1])
        else:
            return (3, 99, item[1])

    final_entries.sort(key=chrom_sort_key)

    header = "track name=IAA41053_181_Input description=InputTargets_AmpliSeqID_IAA41053 solution_type=181 AmpliSeq_Version=8.0 workflow=DNA ionVersion=4.0 db=hg38 reference=GRCh38.p2 cosmicVersion=v86 type=bedDetail color=25,125,255 priority=1\n"

    with open(output_file, 'w', encoding='utf-8') as out:
        out.write(header)
        for chrom, start, end, uid, rsid, score in final_entries:
            out.write(f"{chrom}\t{start}\t{end}\t{uid}\t{rsid}\t{score}\n")

    print("\n" + "=" * 70)
    print("RELATÓRIO COM PROBES DE 120PB - CR-PGX-V9")
    print("=" * 70)
    print(f"Alvos mantidos do seu desenho original: {len(valid_entries)}")
    print(f"Alvos com falha removidos: {len(removed_entries)}")
    print(f"Probes de 120pb adicionadas (Illumina): {len(added_entries)}")
    print(f"Total final no novo BED: {len(final_entries)}")
    print("=" * 70)
    print(f"Arquivo gerado com sucesso: {output_file}\n")


if __name__ == "__main__":
    generate_cr_pgx_v9()