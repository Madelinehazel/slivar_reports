import pandas as pd
import argparse
from parse_report import parse_functions


def main(report):
    report = parse_functions.apply_parse_spliceAI(report)
    report = parse_functions.apply_zygosity(report)
    report = report.drop(columns=["spliceAI_parsed", "spliceai_score"])
    report.to_csv("test.txt", index=False, sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Parses tsv output from slivar into C4R report"
    )
    parser.add_argument("-report", type=str, help="input report text file")

    args = parser.parse_args()

    filename = args.report.strip(".txt")
    report = pd.read_csv(args.report, sep="\t", encoding="utf-8")

    main(report)
