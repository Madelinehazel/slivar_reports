import pandas as pd
import argparse
from parse_report import parse_functions
from datetime import date



def main(report, hpo, filename):
    report = parse_functions.apply_parse_spliceAI(report)
    report = parse_functions.apply_zygosity(report)
    report = parse_functions.apply_add_omim('omim_phenotype', report)
    report = parse_functions.apply_add_omim('omim_inheritance', report)
    #report = parse_functions.apply_alt_depth(report)
    report = parse_functions.apply_ucsc_link(report)
    report = parse_functions.apply_gnomad_link(report)
    report = parse_functions.get_ensemble_gene(report)
    report = parse_functions.add_gnomadscore(report)
    report = parse_functions.add_gene_description(report)
    report = parse_functions.add_orphanet(report)
    report = report.rename({"gene": "Gene"}, axis=1)
    report = parse_functions.add_imprint(report)
    report = parse_functions.add_pseudoautosomal(report)
    report = parse_functions.vest3_score(report)
    report = parse_functions.add_hgmd(report)
    # report = parse_functions.add_c4r_exome_db(report)
    report = parse_functions.apply_parse_enst(report)
    report = parse_functions.apply_parse_consequence(report)
    report = parse_functions.apply_format_highest_impact(report)
    report["HPO"] = [parse_functions.add_hpo(hpo, gene) for gene in report['Gene'].values]

    report = report.rename(
        {
            "gene_impact_transcript_Gene_CANONICAL_EXON_HGVSc_HGVSp_Protein_position_Consequence_PolyPhen_SIFT_DOMAINS": "Info"
        },
        axis=1,
    )
    report = report.rename(
        {"depths(sample,dad,mom)": "trio_coverage(sample,dad,mom)"}, axis=1
    )
    report = parse_functions.merge_variants(report)
    report = report.rename({"clinvar_pathogenic": "Clinvar"}, axis=1)
    report["Clinvar"] = report["Clinvar"].fillna("NA")
    report = report.rename({"gnomad_nhomalt": "gnomad_hom(v2.1.1)"}, axis=1)
    report = report.rename({"phyloP30way_mammalian": "Conserved_in_30_mammals"}, axis=1)
    report["Conserved_in_30_mammals"] = report["Conserved_in_30_mammals"].fillna("NA")
    report = report.rename({"CADD_phred": "Cadd_score"}, axis=1)
    report["Cadd_score"] = report["Cadd_score"].fillna("NA")
    report["REVEL_score"] = report["REVEL_score"].fillna("NA")
    report["Gerp_score"] = report["Gerp_score"].fillna("NA")
    report["UCE_100bp"] = report["UCE_100bp"].fillna("0")
    report["UCE_100bp"] = report["UCE_100bp"].replace("UCE_100bp", "1")
    report["UCE_200bp"] = report["UCE_200bp"].fillna("0")
    report["UCE_200bp"] = report["UCE_100bp"].replace("UCE_200bp", "1")
    for col in ["gnomad_af_popmax", "gnomad_af", "gnomad_ac", "gnomad_hom", "hprc_af", "hprc_ac", "hprc_hom", "cmh_af", "cmh_ac"]:
        report = parse_functions.sub_period_with_zero(report, col)
    report = report.drop(
        columns=[
            "spliceAI_parsed",
            "spliceai_score",
            "PAR_Gene_name",
            "gene_description_1",
            "gene_description_2",
            "external_gene_name",
            "info_parsed",
        ]
    )  # gene_description_1 is pLI_score

    # long read specific filters
    report = parse_functions.filter_hpc_cmh(report)
    report = report[report['chr:pos:ref:alt'].apply(lambda row: parse_functions.filter_linkage_region(row, 'chr3', 178477395 , 183436739))]
    # gnomad_nhomalt from slivar uses gnomadv2.1.1, it's the updated value of gnomad_hom
    # Burden(sample,dad,mom)
    # add Ensembl_transcript_id and info columns
    # Number_of_callers; Old-multiallelic
    report = report[
        [
            "#mode",
            "family_id",
            "sample_id",
            "chr:pos:ref:alt",
            "UCSC_link",
            "GNOMAD_link",
            "zygosity(PB07_CH0076;PB09_CH0627;PB10_CH1073)",
            "Gene",
            "genotype(PB07_CH0076;PB09_CH0627;PB10_CH1073)",
            "Variation",
            "Info",
            "RefSeq_change",
            "RefSeq_canonical",
            #           "DP",
            #           "QUAL",
            # "alt_depth(sample,dad,mom)",
            "trio_coverage(PB07_CH0076;PB09_CH0627;PB10_CH1073)",
            "allele_balance(PB07_CH0076;PB09_CH0627;PB10_CH1073)",
            "Ensembl_gene_id",
            "Gene_description",
            "HPO",
            "omim_phenotype",
            "omim_inheritance",
            "Orphanet",
            "Clinvar",
            # "Frequency_in_C4R",
            # "Seen_in_C4R_samples",
            "HGMD_id",
            "HGMD_gene",
            "HGMD_tag",
            "HGMD_ref",
            "gnomad_af_popmax",
            "gnomad_af",
            "gnomad_ac",
            "gnomad_hom",
            "hprc_af",
            "hprc_ac",
            "hprc_hom",
            "cmh_af",
            "cmh_ac",
            #           "gnomad_hom(v2.1.1)",
            "Ensembl_transcript_id",
            "AA_position",
            "Exon",
            "Protein_domains",
            "rs_ids",
            "Gnomad_oe_lof_score",
            "Gnomad_oe_mis_score",
            "Exac_pli_score",
            "Exac_prec_score",
            "Exac_pnull_score",
            "Conserved_in_30_mammals",
            "spliceAI_impact",
            "spliceAI_score",
            "Polyphen_score",
            "Sift_score",
            "Cadd_score",
            "Vest3_score",
            "REVEL_score",
            "Gerp_score",
            "Imprinting_status",
            "Imprinting_expressed_allele",
            "Pseudoautosomal",
            "UCE_100bp",
            "UCE_200bp",
        ]
    ]

    report = report.fillna('NA')
    today = date.today()
    today = today.strftime("%Y-%m-%d")
    report.to_csv(f"{filename}_parsed_GRCh38_{today}.csv", index=False, sep=",", encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Parses tsv output from slivar into C4R report"
    )
    parser.add_argument("-report", type=str, help="input report text file", required=True)
    parser.add_argument("-hpo", type=str, help="hpo terms")

    args = parser.parse_args()

    filename = args.report.strip(".txt")
    report = pd.read_csv(args.report, sep="\t", encoding="utf-8")
    hpo = pd.read_csv(args.hpo, comment="#", skip_blank_lines=True, sep='\t', encoding="ISO-8859-1", engine="python")
    hpo.columns = hpo.columns.str.strip()


    main(report, hpo, filename)
