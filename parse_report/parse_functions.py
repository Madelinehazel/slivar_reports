import pandas as pd
import os
import re

default_tables_path = "~/cre/data/"
hgmd_filepath = "/hpf/largeprojects/ccm_dccforge/dccforge/results/database/hgmd.csv"


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


# OMIM inheritance
# OMIM phenotype
OMIM_file = os.path.join(
    default_tables_path, "OMIM_hgnc_join_omim_phenos_2021-04-27.tsv"
)
OMIM = pd.read_csv(OMIM_file, sep="\t")
OMIM = OMIM.rename(columns={"gene_name": "gene"})
OMIM = OMIM.dropna(subset=['omim_phenotype'])


def add_omim(col, gene):
    if pd.isna(gene):
        return "NA"
    else:
        omim = []
        gene = gene.split(";")
        for g in gene:
            try:
                omim.append(
                    str(OMIM[OMIM["gene"] == g][col].values[0])
                )
            except IndexError:
                pass
        omim = ",".join(omim)
        return omim


def apply_add_omim(col, report):
    report[col] = report["gene"].apply(lambda x: add_omim(col, x))
    report[col] = report[col].fillna("NA")
    report[col] = report[col].replace("", "NA")
    return report


def alt_depth(DP, AB):
    DP = DP.split(",")
    DP = [float(num) for num in DP]
    AB = AB.split(",")
    AB = [float(num) for num in AB]
    AD = [a * b for a, b in zip(DP, AB)]
    AD = [int(num) for num in AD]
    AD = ",".join([str(num) for num in AD])

    return AD


def apply_alt_depth(report):
    report["alt_depth(sample,dad,mom)"] = ""

    alt_dp = []
    for index, row in report.iterrows():
        dp = row["depths(sample,dad,mom)"]
        ab = row["allele_balance(sample,dad,mom)"]
        alt_dp.append(alt_depth(dp, ab))
        
    report["alt_depth(sample,dad,mom)"] = alt_dp

    return report


def gnomad_link(info):
    """Generate website link with the specified format. for eg. "1:1581136:G:A" """
    info = info.split(":")
    chr = info[0]
    pos = info[1]
    ref = info[2]
    alt = info[3]
    return '=HYPERLINK("http://gnomad.broadinstitute.org/variant/{chr}-{pos}-{ref}-{alt}","GNOMAD_link")'.format(
        chr=chr, pos=pos, ref=ref, alt=alt
    )


def apply_gnomad_link(report):
    report["GNOMAD_link"] = report["chr:pos:ref:alt"].apply(
        lambda row: gnomad_link(row)
    )

    return report


def ucsc_link(info):
    info = info.split(":")
    chr = info[0]
    pos = info[1]
    ref = info[2]
    alt = info[3]
    return '=HYPERLINK("http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&hgt.out3=10x&position={chr}:{pos}","UCSC_link")'.format(
        chr=chr, pos=pos
    )


def apply_ucsc_link(report):
    info = report["chr:pos:ref:alt"]
    report["UCSC_link"] = [
        ucsc_link(report["chr:pos:ref:alt"][row]) for row in range(len(info))
    ]

    return report


def get_ensemble_gene(report):
    consequences = report[
        "gene_impact_transcript_Gene_CANONICAL_EXON_HGVSc_HGVSp_Protein_position_Consequence_PolyPhen_SIFT_DOMAINS"
    ]
    # ensg = re.findall("ENSG[0-9]+", consequences)[0]
    report["Ensembl_gene_id"] = [
        'None' if re.findall("ENSG[0-9]+", consequences[row])==[] else re.findall("ENSG[0-9]+", consequences[row])[0]
        for row in range(len(consequences))
    ]

    return report


# Gnomad_oe_lof_score
# Gnomad_oe_mis_score
gnomad_score_file = os.path.join(default_tables_path, "gnomad_scores.csv")
gnomad_score = pd.read_csv(gnomad_score_file, sep=",")
gnomad_score.index = list(gnomad_score.index)  # change rangeindex to int64index


def add_gnomadscore(report):
    report = pd.merge(report, gnomad_score, how="left")
    report["Gnomad_oe_lof_score"] = report["Gnomad_oe_lof_score"].fillna("NA")
    report["Gnomad_oe_mis_score"] = report["Gnomad_oe_mis_score"].fillna("NA")
    report["Exac_pli_score"] = report["Exac_pli_score"].fillna("NA")
    report["Exac_prec_score"] = report["Exac_prec_score"].fillna("NA")
    report["Exac_pnull_score"] = report["Exac_pnull_score"].fillna("NA")
    return report


# Gene description
gene_description_file = os.path.join(default_tables_path, "ensembl_w_description.txt")
gene_description = pd.read_csv(gene_description_file, sep="\t")
gene_description = gene_description.rename(
    columns={"ensembl_gene_id": "Ensembl_gene_id"}
)


def add_gene_description(report):
    report = pd.merge(report, gene_description, how="left")
    return report


# orphanet
orphanet_file = os.path.join(default_tables_path, "orphanet.txt")
orphanet = pd.read_csv(orphanet_file, sep="\t", encoding="latin1")


def add_orphanet(report):
    report = pd.merge(report, orphanet, how="left")
    report["Orphanet"] = report["Orphanet"].fillna("NA")
    return report


# Imprinting_status
# Imprinting_expressed_allele
imprinting_file = os.path.join(default_tables_path, "imprinting.txt")
imprinting = pd.read_csv(imprinting_file, sep="\t")


def add_imprint(report):
    report = pd.merge(report, imprinting, how="left")
    report["Imprinting_status"] = report["Imprinting_status"].fillna("NA")
    report["Imprinting_expressed_allele"] = report[
        "Imprinting_expressed_allele"
    ].fillna("NA")
    return report


# pseudoautosomal
pseudoautosomal_file = os.path.join(default_tables_path, "pseudoautosomal.txt")
pseudoautosomal = pd.read_csv(pseudoautosomal_file, sep="\t")


def add_pseudoautosomal(report):
    report = pd.merge(report, pseudoautosomal, how="left")
    report["Pseudoautosomal"] = report["Pseudoautosomal"].fillna("NA")
    return report




# Number of callers
### add a column (placeholder) called "Number_of_callers" and put "Number_of_calleres"


# Burden


# old multiallelic


# HGMD annotation
### add placeholder for HGMD_id, _gene, _tag, _ref with "NA"
hgmd = pd.read_csv(hgmd_filepath, sep=",", header=None, dtype=str)
hgmd.columns = [
    "chrom",
    "pos",
    "HGMD_id",
    "ref",
    "alt",
    "HGMD_gene",
    "HGMD_tag",
    "author",
    "allname",
    "vol",
    "page",
    "year",
    "pmid",
    "rsid",
]
superindex = hgmd["chrom"] + ":" + hgmd["pos"] + ":" + hgmd["ref"] + ":" + hgmd["alt"]
ref = (
    hgmd["author"]
    + ", "
    + hgmd["allname"]
    + ", vol"
    + hgmd["vol"]
    + ", page"
    + hgmd["page"]
    + ","
    + hgmd["year"]
    + ","
    + "PMID:"
    + hgmd["pmid"]
)
refdf = pd.DataFrame(
    {
        "chr:pos:ref:alt": superindex,
        "HGMD_id": hgmd["HGMD_id"],
        "HGMD_gene": hgmd["HGMD_gene"],
        "HGMD_tag": hgmd["HGMD_tag"],
        "HGMD_ref": ref,
    }
)


def add_hgmd(report):
    report = pd.merge(report, refdf, how="left")
    report["HGMD_id"] = report["HGMD_id"].fillna("NA")
    report["HGMD_gene"] = report["HGMD_gene"].fillna("NA")
    report["HGMD_tag"] = report["HGMD_tag"].fillna("NA")
    report["HGMD_ref"] = report["HGMD_ref"].fillna("NA")
    return report


# frequency_in_C4R, seen_in_C4R_samples
def large_sample(sample): # sample is a string in the list
    if len(sample.split("; ")) >= 100:
        first100 = sample.split("; ")[0:99]
        sample = f'{"; ".join(item for item in first100)} ...'
    else:
        pass
    return sample


def add_c4r_exome_db(report):
    c4r_counts_file = "/hpf/largeprojects/ccm_dccforge/dccforge/results/database/seen_in_c4r_counts.txt"
    c4r_sample_file = "/hpf/largeprojects/ccm_dccforge/dccforge/results/database/seen_in_c4r_samples.txt"
    c4r_counts = pd.read_csv(
        c4r_counts_file,
        sep="\t",
        header=0,
        names=["chr:pos:ref:alt", "Frequency_in_C4R"],
        skiprows=0,
        encoding="utf-8",
    )
    c4r_sample = pd.read_csv(
        c4r_sample_file,
        sep="\t",
        header=0,
        names=["chr:pos:ref:alt", "Seen_in_C4R_samples"],
        skiprows=0,
        encoding="utf-8",
    )
    c4r_counts["chr:pos:ref:alt"] = c4r_counts["chr:pos:ref:alt"].str.replace("-", ":")
    c4r_sample["chr:pos:ref:alt"] = c4r_sample["chr:pos:ref:alt"].str.replace("-", ":")
    report = pd.merge(report, c4r_counts, on="chr:pos:ref:alt", how="left")
    report = pd.merge(report, c4r_sample, on="chr:pos:ref:alt", how="left")
    report["Frequency_in_C4R"] = report["Frequency_in_C4R"].fillna("0")
    report["Seen_in_C4R_samples"] = report["Seen_in_C4R_samples"].fillna("0")

    #sample_list = report["Seen_in_C4R_samples"].tolist()
    report["Seen_in_C4R_samples"] = report["Seen_in_C4R_samples"].apply(
        lambda row: large_sample(row)
    )
    return report

# parse Ensembl transcript
def parse_enst(consequence):
    consequence = consequence.split(";")
    enst = [transcript for transcript in consequence if re.findall("ENST[0-9]+", transcript.split("/")[2])]
    if len(enst) == 0:
        result = ["NA"]

    else:
        enst_yes =  [transcript for transcript in enst if re.match("YES", transcript.split("/")[4])]
        if len(enst_yes) == 0:
            result = enst[0].split("/")[2]
        elif len(enst_yes) >= 1:
            result = enst_yes[0].split("/")[2]
    
    return result

def apply_parse_enst(report):
    report["Ensembl_transcript_id"] = report['gene_impact_transcript_Gene_CANONICAL_EXON_HGVSc_HGVSp_Protein_position_Consequence_PolyPhen_SIFT_DOMAINS'].apply(
        lambda row: parse_enst(row)
    )
    return report

# parse refseq
def get_anno_list(anno): # after split("/"), used for parse_inf() function
    #gene = anno[0]
    #nm = anno[2]
    #when there is HGVC coding and protein information:
    #canonical = anno[4]
    #exon = anno[5:6]
    #refseq_HGVC = anno[7]
    #refseq_HGVCp = anno[8]
    #protein = anno[9:10]
    #protein_domain = anno[14]
    #polyphen = anno[12]
    #sift = anno[13]
    
    exon = anno[5]
    if len(exon) == 0:
        result = [anno[0], anno[2], anno[6], anno[7], exon, "NA", anno[12], anno[10], anno[11]]
                
            
    else:
        exon = f"{anno[5]}_{anno[6]}"
        if anno[9] == "":
            result = [anno[0], anno[2], anno[7],  anno[8], exon, "NA", anno[13], anno[11], anno[12]]
        else:
            protein = f"{anno[9]}_{anno[10]}"
            result = [anno[0], anno[2], anno[7],  anno[8], exon, protein, anno[14], anno[12], anno[13]]
    return result

def parse_inf(consequence):
    consequence = consequence.split(";")
    refseq = [transcript for transcript in consequence if re.findall("NM_[0-9]+", transcript.split("/")[2])]
    
    if len(refseq) == 0:
        result = ["NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"]
        # result = [refseq_HGVC, HGVCp, exon, protein, protein_domain, polyphen, sift, canonical]        
            
    else:
        nm_yes = [transcript for transcript in refseq if re.match("YES", transcript.split("/")[4])]
        if len(nm_yes) == 0:   
            canonical = "No"
            anno = refseq[0].split("/")    # take the first refseq transcript when there is no canonical
            result = get_anno_list(anno)
            result.append(canonical)
            result = list(map(get_change, result))
            coding = ":".join(result[1:3])
            np = result[3]
            result[0] = ":".join([result[0], coding, np]) # gene:NM_:coding_change:protein_change
                
              

        else:    # there is NM_ and YES
            canonical = "Yes"
            if len(nm_yes) == 1: # there is exactly one NM_ canonical
                anno = nm_yes[0].split("/")
                result = get_anno_list(anno)
                result.append(canonical)
                result = list(map(get_change, result))
                coding = ":".join(result[1:3])
                np = result[3]
                result[0] = ":".join([result[0], coding, np]) # gene:NM_:coding_change:protein_change

                
            elif len(nm_yes) > 1: # there are > 1 NM_ canonical, could be in different genes
                parseresult=[]
                for refcano in nm_yes:
                    anno = refcano.split("/")
                    parseresult.append(get_anno_list(anno)) # list of lists without cano yet
                
                parseresult_refseq = format_refseq(parseresult)
                result = list((parseresult_refseq, parseresult[0][1:], canonical)) # use the first refseq cano to extract all other information
                result = list(iterFlatten(result))
                result[1:] = list(map(get_change, result[1:]))             

    return result
    
def apply_parse_consequence(report):
    report["info_parsed"] = report['gene_impact_transcript_Gene_CANONICAL_EXON_HGVSc_HGVSp_Protein_position_Consequence_PolyPhen_SIFT_DOMAINS'].apply(
        lambda row: parse_inf(row)
    )
    refseq_change = []
    nm = []
    refseq_hgvc = []
    np = []
    exon = []
    protein = []
    protein_domain = []
    polyphen = []
    sift = []
    canonical = []

    for index, row in report.iterrows():
        parsed = row['info_parsed'] # first row
       
        refseq_change.append(parsed[0])
        nm.append(parsed[1])
        refseq_hgvc.append(parsed[2]) # first element in the list
        np.append(parsed[3])
        exon.append(parsed[4])
        protein.append(parsed[5])
        protein_domain.append(parsed[6])
        polyphen.append(parsed[7])
        sift.append(parsed[8])
        canonical.append(parsed[9])
    
    
    report["RefSeq_change"] = refseq_change
    #report["nm"] = nm
    #report["refseq_hgvc"] = refseq_hgvc
    #report["np"] = np
    report["Exon"] = exon
    report["AA_position"] = protein
    report["Protein_domains"] = protein_domain
    report["Polyphen_score"] = polyphen
    report["Sift_score"] = sift
    report["RefSeq_canonical"] = canonical

    report["Protein_domains"] = report["Protein_domains"].replace("", "NA")
    report["Polyphen_score"] = report["Polyphen_score"].replace("", "None")
    report["Sift_score"] = report["Sift_score"].replace("", "None")
    
    return report

### sub functions for multiple refseq_cano transcripts
def format_refseq(parseresult):
    refseq = []
    for transcript in parseresult:
        coding = get_hgvc_coding(transcript)
        protein_change = get_hgvc_protein(transcript)
        refseq.append([transcript[0], coding, protein_change])
    
    result = []
    for trans in range(len(refseq)): 
        result.append(
        ":".join(refseq[trans])
        )    
    result = ",".join(result)

    return result

def get_hgvc_coding(annotranscript):
    if annotranscript[2] == "":
        return ":".join([annotranscript[1], "NA"]) # has NM_ only
    else:
        return annotranscript[2]

def get_hgvc_protein(annotranscript):
    if annotranscript[3] == "":
        return "NA"
    else:
        return annotranscript[3].split(":")[1]

def iterFlatten(root):    # function to flatten lists of list
    if isinstance(root, (list, tuple)):
        for element in root:
            for e in iterFlatten(element):
                yield e
    else:
        yield root

##################################
def get_change(change):
    if change == "":
        result_change = "NA"
    elif change == "NA":
        result_change = "NA"
    elif "NP_" in change:
        result_change = change.split(":")[1]
    elif "NM_" and ":c" in change:
        result_change = change.split(":")[1]    
    else:
        result_change = change
    return result_change



### FIX FORMATING IN THE REPORT

# Column1 - Mode of inheritance
# Column2 - FamilyID
# Column3 - SampleID
# Column4 - chr:pos:ref:alt
# Column5 - genotype(sample,dad,mom) (must be placed after zygosity function)

# Column6 - DP
# Column7 - clinvar_pathogenic

# Column8 - gnmoad_af_popmax (from vcfanno)

# Column10 - gnomad_af
# Column11 - gnomad_ac
# Column12 - gnomad_hom
# Column13 - gnomad_nhomalt
def gnomad_replace(report):
    report["gnomad_af_popmax"] = report["gnomad_af_popmax"].replace(".", 0)
    report["gnomad_af"] = report["gnomad_af"].replace(".", 0)
    report["gnomad_ac"] = report["gnomad_ac"].replace(".", 0)
    report["gnomad_hom"] = report["gnomad_hom"].replace(".", 0)
    return report


# Vest3_score
def vest3_score(report):
    report["Vest3_score"] = report["Vest3_score"].apply(
        lambda row: max((str.split(str(row), sep=",")))
    )
    report["Vest3_score"] = report["Vest3_score"].replace("nan", "None")
    return report

# format highes_impact column, renamed to "Variation"
def format_highest_impact(impact):
    rank_impact = impact.split("_")
    impact = " ".join(rank_impact[1:])
    return impact

def apply_format_highest_impact(report):
    report["Variation"] = report["highest_impact"].apply(
        lambda row: format_highest_impact(row)
    )
    return report

# module load python/3.7.1
# python3 generate_report.py -report ../report/report.txt
