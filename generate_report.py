import pandas as pd
import argparse
from parse_report import parse_functions
import os


def main(report):
    report = parse_functions.apply_parse_spliceAI(report)
    report = parse_functions.apply_zygosity(report)
    report = parse_functions.add_omim(report)
    report = parse_functions.apply_alt_depth(report)
    report = parse_functions.apply_ucsc_link(report)
    report = parse_functions.apply_gnomad_link(report)
    report = parse_functions.get_ensemble_gene(report)
    report = parse_functions.add_gnomadscore(report)
    report = parse_functions.add_gene_description(report)
    report = parse_functions.add_orphanet(report)
    report = report.rename({'gene': 'Gene'}, axis=1)
    report = parse_functions.add_imprint(report)
    report = parse_functions.add_pseudoautosomal(report)
    report = report.drop(columns=["spliceAI_parsed", "spliceai_score", "PAR_Gene_name"])
    report = report.rename({'gene_description_1': 'pLI_score'}, axis=1)
    report = report.rename({'gene_description_2': 'gene_description'}, axis=1)
    report = report.rename({'gene_impact_transcript_Gene_EXON_HGVSc_HGVSp_Protein_position_Consequence_PolyPhen_SIFT': 'Info'}, axis=1)
    report.to_csv("test.txt", index=False, sep="\t", quotechar="'", escapechar="\'", doublequote=False)
    report.to_csv("test.csv", index=False, sep=",", quotechar="'", escapechar="\'", doublequote=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Parses tsv output from slivar into C4R report"
    )
    parser.add_argument("-report", type=str, help="input report text file")

    args = parser.parse_args()

    filename = args.report.strip(".txt")
    report = pd.read_csv(args.report, sep="\t", encoding="utf-8")

    

    main(report)
