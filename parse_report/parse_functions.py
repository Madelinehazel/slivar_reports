import pandas as pd
import numpy as np
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
OMIM_file = os.path.join(default_tables_path, "OMIM_hgnc_join_omim_phenos_2021-04-27.tsv")
OMIM = pd.read_csv(OMIM_file, sep = "\t")
OMIM = OMIM.rename(columns={"gene_name": "gene"})
def add_omim(report):
    report = pd.merge(report, OMIM, how = "left")
    report['omim_phenotype'] = report['omim_phenotype'].fillna('NA')
    report['omim_inheritance'] = report['omim_inheritance'].fillna('NA')
    return report


def alt_depth(DP, AB):
    DP = DP.str.split(pat=",")
    DP = np.asarray(DP.tolist(), dtype='float32')
    AB = AB.str.split(pat=",")
    AB = np.asarray(AB.tolist(), dtype='float32')
    AD=np.array([DP[:,0] * AB[:,0], DP[:,1] * AB[:,1], DP[:,2] * AB[:,2]])
    AD = AD.transpose().astype(int)
    AD = AD.astype('str')
    
    return AD

def apply_alt_depth(report):
    DP = report["depths(sample,dad,mom)"]
    AB = report["allele_balance(sample,dad,mom)"]
    
    AD_column = alt_depth(DP, AB)
    
    report["alt_depth(sample,dad,mom)"] = [','.join(row) for row in AD_column]
    
    return report

def gnomad_link(info):
    """Generate website link with the specified format. for eg. "1:1581136:G:A"  """
    info=info.split(':')
    chr=info[0]
    pos=info[1]
    ref=info[2]
    alt=info[3]
    return "=HYPERLINK(\"http://gnomad.broadinstitute.org/variant/{chr}-{pos}-{ref}-{alt}\",\"GNOMAD_link\")".format(chr=chr, pos=pos, ref=ref, alt=alt)

def apply_gnomad_link(report):
    report["GNOMAD_link"] = report["chr:pos:ref:alt"].apply(
        lambda row: gnomad_link(row)
    )

    return report


def ucsc_link(info):
    info=info.split(':')
    chr=info[0]
    pos=info[1]
    ref=info[2]
    alt=info[3]
    return "=HYPERLINK(\"http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&hgt.out3=10x&position={chr}:{pos}\",\"UCSC_link\")".format(chr=chr, pos=pos)

def apply_ucsc_link(report):
    info=report["chr:pos:ref:alt"]
    report["UCSC_link"] = [ucsc_link(report["chr:pos:ref:alt"][row]) for row in range(len(info))]

    return report



def get_ensemble_gene(report):
    consequences = report['gene_impact_transcript_Gene_CANONICAL_EXON_HGVSc_HGVSp_Protein_position_Consequence_PolyPhen_SIFT_DOMAINS']
    #ensg = re.findall("ENSG[0-9]+", consequences)[0]
    report["Ensembl_gene_id"] = [re.findall("ENSG[0-9]+", consequences[row])[0] for row in range(len(consequences))]

    return report

# Gnomad_oe_lof_score
# Gnomad_oe_mis_score
gnomad_score_file = os.path.join(default_tables_path, "gnomad_scores.csv")
gnomad_score = pd.read_csv(gnomad_score_file, sep=",")
gnomad_score.index=list(gnomad_score.index) # change rangeindex to int64index
def add_gnomadscore(report):
    report = pd.merge(report, gnomad_score, how="left")
    report['Gnomad_oe_lof_score'] = report['Gnomad_oe_lof_score'].fillna('NA')
    report['Gnomad_oe_mis_score'] = report['Gnomad_oe_mis_score'].fillna('NA')
    report['Exac_pli_score'] = report['Exac_pli_score'].fillna('NA')
    report['Exac_prec_score'] = report['Exac_prec_score'].fillna('NA')
    report['Exac_pnull_score'] = report['Exac_pnull_score'].fillna('NA')
    return report


# Gene description
gene_description_file = os.path.join(default_tables_path, "ensembl_w_description.txt")
gene_description = pd.read_csv(gene_description_file, sep="\t")
gene_description = gene_description.rename(columns={"ensembl_gene_id": "Ensembl_gene_id"})
def add_gene_description(report):
    report = pd.merge(report, gene_description, how = "left")
    return report


# orphanet 
orphanet_file = os.path.join(default_tables_path, "orphanet.txt")
cols=pd.read_csv(orphanet_file, sep="\t").columns
orphanet= pd.read_csv(orphanet_file, sep="\t", usecols=cols[0:2])
def add_orphanet(report):
    report = pd.merge(report, orphanet, how = "left")
    report['Orphanet'] = report['Orphanet'].fillna('NA')
    return report



# Imprinting_status
# Imprinting_expressed_allele
imprinting_file = os.path.join(default_tables_path, "imprinting.txt")
imprinting = pd.read_csv(imprinting_file, sep = "\t")
def add_imprint(report):
    report = pd.merge(report, imprinting, how = "left")
    report['Imprinting_status'] = report['Imprinting_status'].fillna('NA')
    report['Imprinting_expressed_allele'] = report['Imprinting_expressed_allele'].fillna('NA')
    return report

# pseudoautosomal
pseudoautosomal_file = os.path.join(default_tables_path, "pseudoautosomal.txt")
pseudoautosomal = pd.read_csv(pseudoautosomal_file, sep = "\t")
def add_pseudoautosomal(report):
    report = pd.merge(report, pseudoautosomal, how = "left")
    report['Pseudoautosomal'] = report['Pseudoautosomal'].fillna('NA')
    return report

# Refseq_change
### add a column (placeholder) called "Refseq_change" and put "NA" in 
### there is no ref_seq notation in the info generated by slivar
def add_placeholder(report, column_name, placeholder):
    report[column_name] = placeholder
    return report


# Number of callers
### add a column (placeholder) called "Number_of_callers" and put "Number_of_calleres"


# Burden


# old multiallelic



# HGMD annotation
### add placeholder for HGMD_id, _gene, _tag, _ref with "NA"
hgmd = pd.read_csv(hgmd_filepath, sep = ",", header=None, dtype=str)
hgmd.columns = ["chrom","pos","HGMD_id","ref","alt","HGMD_gene","HGMD_tag","author",
                          "allname","vol","page","year","pmid","rsid"]
superindex = hgmd["chrom"] + ":" + hgmd["pos"] + ":" + hgmd["ref"] + ":" + hgmd["alt"]
ref = hgmd['author'] + ', ' + hgmd['allname'] + ', vol' + hgmd['vol'] + ', page' + hgmd['page'] + ',' + hgmd['year'] + ',' + 'PMID:' + hgmd['pmid']
refdf = pd.DataFrame({'chr:pos:ref:alt': superindex, 'HGMD_id': hgmd['HGMD_id'], 'HGMD_gene': hgmd['HGMD_gene'],  'HGMD_tag': hgmd['HGMD_tag'], 'HGMD_ref': ref})


def add_hgmd(report):
    report = pd.merge(report, refdf, how = "left")
    report['HGMD_id'] = report['HGMD_id'].fillna('NA')
    report['HGMD_gene'] = report['HGMD_gene'].fillna('NA')
    report['HGMD_tag'] = report['HGMD_tag'].fillna('NA')
    report['HGMD_ref'] = report['HGMD_ref'].fillna('NA')
    return report


# frequency_in_C4R, seen_in_C4R_samples


### FIX FORMATING IN THE REPORT
    
#Column1 - Mode of inheritance
#Column2 - FamilyID
#Column3 - SampleID
#Column4 - chr:pos:ref:alt
#Column5 - genotype(sample,dad,mom) (must be placed after zygosity function)

#Column6 - DP
#Column7 - clinvar_pathogenic

#Column8 - gnmoad_af_popmax (from vcfanno)

#Column10 - gnomad_af
#Column11 - gnomad_ac
#Column12 - gnomad_hom
#Column13 - gnomad_nhomalt
def gnomad_replace(report):
    report['gnomad_af_popmax'] = report['gnomad_af_popmax'].replace('.', 0)
    report['gnomad_af'] = report['gnomad_af'].replace('.', 0)
    report['gnomad_ac'] = report['gnomad_ac'].replace('.', 0)
    report['gnomad_hom'] = report['gnomad_hom'].replace('.', 0)
    return report

# Vest3_score
def vest3_score(report):
    report['Vest3_score'] = report['Vest3_score'].apply(
        lambda row: max((str.split(str(row), sep = ',')))
        )
    report['Vest3_score'] = report['Vest3_score'].replace('nan', 'None')
    return report





# python3 generate_report.py -report ../report/report.txt