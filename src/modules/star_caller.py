import os
import gzip
import pysam
from cyvcf2 import VCF

def call_star_alleles_batch_from_vcf(vcf_path: str, gene_symbols: list[str] = None) -> list[dict]:
    # Verifica se já existe a versão ordenada e indexada para evitar reprocessamento desnecessário
    sorted_vcf_path = vcf_path.replace(".vcf.gz", ".sorted.vcf.gz")
    if not sorted_vcf_path.endswith('.vcf.gz'):
        sorted_vcf_path += ".vcf.gz"

    tbi_path = sorted_vcf_path + ".tbi"

    if not os.path.exists(tbi_path):
        print("VCF precisa de normalização e ordenação. Processando via Python...")
        header = []
        variants_data = []

        with gzip.open(vcf_path, 'rt', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if line.startswith('#'):
                    header.append(line)
                else:
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        chrom = parts[0]
                        try:
                            pos = int(parts[1])
                            variants_data.append((chrom, pos, line))
                        except ValueError:
                            continue

        chrom_order = {f"chr{i}": i for i in range(1, 23)}
        chrom_order.update({str(i): i for i in range(1, 23)})
        chrom_order.update({'chrX': 23, 'X': 23, 'chrY': 24, 'Y': 24, 'chrMT': 25, 'MT': 25, 'M': 25})

        variants_data.sort(key=lambda x: (chrom_order.get(x[0], 99), x[1]))

        temp_sorted_txt = vcf_path.replace(".vcf.gz", ".sorted.vcf")
        with open(temp_sorted_txt, 'w', encoding='utf-8') as out:
            for h in header:
                out.write(h)
            for _, _, line in variants_data:
                out.write(line if line.endswith('\n') else line + '\n')

        pysam.tabix_compress(temp_sorted_txt, sorted_vcf_path, force=True)
        pysam.tabix_index(sorted_vcf_path, preset="vcf", force=True)

        if os.path.exists(temp_sorted_txt):
            os.remove(temp_sorted_txt)

    vcf_target = sorted_vcf_path if os.path.exists(sorted_vcf_path) else vcf_path

    queryset = Gene.objects.all()
    if gene_symbols:
        queryset = queryset.filter(name__in=[g.upper() for g in gene_symbols])

    queryset = queryset.prefetch_related('haplotypes__variants', 'variant_set')

    vcf = VCF(vcf_target)
    results = []

    for gene in queryset:
        chrom = gene.chromosome
        all_gene_variants = list(gene.variant_set.all())
        if not all_gene_variants:
            print(f"[{gene.name}] Aviso: Nenhuma variante cadastrada na tabela Variant para este gene no banco.")
            continue

        start = min(int(v.position) for v in all_gene_variants if v.position)
        end = max(int(v.position) for v in all_gene_variants if v.position)

        variants_in_patient = set()
        try:
            query_str = f"{chrom}:{start}-{end}"
            for variant in vcf(query_str):
                for alt in variant.ALT:
                    key = (variant.POS, variant.REF, alt)
                    raw_genos = variant.genotypes
                    if raw_genos and len(raw_genos) > 0:
                        sample_alleles = raw_genos[0][:2]
                        if any(g is not None and g > 0 for g in sample_alleles):
                            variants_in_patient.add(key)
        except Exception as e:
            print(f"Erro na query {query_str}: {e}")

        print(f"\n--- Diagnóstico para {gene.name} ({chrom}:{start}-{end}) ---")
        print(f"Variantes mutadas encontradas no VCF do paciente: {list(variants_in_patient)}")

        matched_haplotypes = []
        haplotypes = sorted(gene.haplotypes.all(), key=lambda h: 0 if h.name != '*1' else 1)

        for hap in haplotypes:
            hap_variants = hap.variants.all()
            if not hap_variants.exists():
                continue

            match = True
            missing_vars = []
            for hv in hap_variants:
                if not hv.position:
                    match = False
                    break
                hv_key = (hv.position, hv.ref, hv.alt)
                if hv_key not in variants_in_patient:
                    match = False
                    missing_vars.append(hv_key)
                    break

            if match:
                matched_haplotypes.append(hap.name)
            else:
                print(f"  Haplótipo {hap.name} descartado. Faltou a variante: {missing_vars}")

        print(f"Haplótipos compatíveis detectados: {matched_haplotypes}")

        if len(matched_haplotypes) == 0:
            matched_hap1 = "*1"
            matched_hap2 = "*1"
        elif len(matched_haplotypes) == 1:
            matched_hap1 = matched_haplotypes[0]
            matched_hap2 = "*1"
        else:
            matched_hap1 = matched_haplotypes[0]
            matched_hap2 = matched_haplotypes[1]

        activity_score = 2.0
        if "*1" not in [matched_hap1, matched_hap2]:
            activity_score = 1.0

        results.append({
            "gene": gene.name,
            "chromosome": gene.chromosome,
            "hap1": matched_hap1,
            "hap2": matched_hap2,
            "diplotype": f"{matched_hap1}/{matched_hap2}",
            "activity_score": activity_score
        })

    vcf.close()
    return results


if __name__ == '__main__':
    vcf_file = "/home/esther/GitProjects/flexlab/backend/genomics/tests/f32_test1/file1.vcf.gz"
    genes_interesse = ["CYP2D6", "CYP2C19"]

    resultados = call_star_alleles_batch_from_vcf(vcf_file, genes_interesse)
    print("\nResultados Finais:")
    for res in resultados:
        print(res)