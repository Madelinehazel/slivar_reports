import pandas as pd


def zygosity(genotype):
    """converts numeric genotype to zygosity for ease of interpretation"""
    zygosity_dict = {"0": "-", "1": "het", "2": "hom", "-1": "missing"}
    zygosity = [zygosity_dict[gt] for gt in genotype.split(",")]
    return (",").join(zygosity)


def apply_zygosity(report):
    report["zygosity(sample,dad,mom)"] = report["genotype(sample,dad,mom)"].apply(
        lambda row: zygosity(row)
    )
    return report


def parse_spliceAI(spliceAI):
    """extracts maximum splice score and associated impact from spliceAI prediction for ease of filtering
    inputs follow the format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL
    example input: T|RYR1|0.00|0.00|0.91|0.08|-28|-46|-2|-31
    example output:  ["RYR1|donor_gain|-2", 0.91]"""
    if pd.isna(spliceAI):
        return ["NA|NA|NA", 0]
    else:
        spliceai = spliceAI.split(",")
        score_dict = {"gene": "NA", "impact": "NA", "score": 0, "pos": "NA"}
        for anno in spliceai:
            # can be more than one spliceAI annotation per variant if it falls in more than one gene
            anno = anno.split("|")
            gene = anno[1]
            DS_AG = anno[2]
            DS_AL = anno[3]
            DS_DG = anno[4]
            DS_DL = anno[5]
            DP_AG = anno[6]
            DP_AL = anno[7]
            DP_DG = anno[8]
            DP_DL = anno[9]
            scores = {
                "acceptor_gain": float(DS_AG),
                "acceptor_loss": float(DS_AL),
                "donor_gain": float(DS_DG),
                "donor_loss": float(DS_DL),
            }
            # find the impact with the maximum delta score
            max_score = 0
            for impact, score in scores.items():
                if score > max_score:
                    max_score = score
            max_score_impact = max(scores, key=scores.get)
            # if all scores are 0, there is no splice impact
            if max_score_impact == 0:
                impact = "NA"
            if score_dict["score"] < max_score:
                score_dict["score"] = max_score
                score_dict["gene"] = gene
                score_dict["impact"] = max_score_impact
                # get the position of the splice change with max predicted delta score
                if max_score_impact == "acceptor_gain":
                    score_dict["pos"] = DP_AG
                elif max_score_impact == "acceptor_loss":
                    score_dict["pos"] = DP_AL
                elif max_score_impact == "donor_gain":
                    score_dict["pos"] = DP_DG
                else:
                    score_dict["pos"] = DP_DL
        return [
            f"{score_dict['gene']}|{score_dict['impact']}|{score_dict['pos']}",
            score_dict["score"],
        ]


def apply_parse_spliceAI(report):
    """splits the impact and score generate by parse_spliceAI into two columns
    not very elegant, but I want to avoid passing the whole df row to parse_spliceAI so
    it is easier to test"""
    report["spliceAI_parsed"] = report["spliceai_score"].apply(
        lambda row: parse_spliceAI(row)
    )
    impact = []
    score = []
    for index, row in report.iterrows():
        parsed = row["spliceAI_parsed"]
        score.append(parsed[1])
        impact.append(parsed[0])
    report["spliceAI_impact"] = impact
    report["spliceAI_score"] = score
    return report
