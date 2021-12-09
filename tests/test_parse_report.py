import pytest
from parse_report import parse_functions


def test_parse_spliceAI():
    AG = parse_functions.parse_spliceAI("T|PGM1|0.23|0.13|0.00|0.00|3|-33|-33|28")
    AL = parse_functions.parse_spliceAI("T|MORN1|0.00|0.08|0.00|0.00|-13|44|44|-23")
    DG = parse_functions.parse_spliceAI("T|RYR1|0.00|0.00|0.91|0.08|-28|-46|-2|-31")
    DL = parse_functions.parse_spliceAI("G|NIT1|0.00|0.00|0.08|0.27|1|6|32|1")
    zero = parse_functions.parse_spliceAI("C|MAST2|0.00|0.00|0.00|0.00|4|-1|19|-1")
    multiple = parse_functions.parse_spliceAI(
        "C|KTI12|0.00|0.00|0.04|0.00|1|3|1|-3,C|TXNDC12|0.00|0.00|0.00|0.00|1|-3|1|-3"
    )

    assert AG == ["PGM1|acceptor_gain|3", 0.23]
    assert AL == ["MORN1|acceptor_loss|44", 0.08]
    assert DG == ["RYR1|donor_gain|-2", 0.91]
    assert DL == ["NIT1|donor_loss|1", 0.27]
    assert zero == ["NA|NA|NA", 0]
    assert multiple == ["KTI12|donor_gain|1", 0.04]


def test_zygosity():
    assert parse_functions.zygosity("1,0,0") == "het,-,-"
    assert parse_functions.zygosity("2,1,0") == "hom,het,-"
    assert parse_functions.zygosity("2,1,-1") == "hom,het,missing"


def test_alt_depth():
    # DP="84,72,69"
    # AB="0.238095,0,0.246377"
    assert parse_functions.alt_depth("84,72,69", "0.238095,0,0.246377") == '19,0,17'
    assert parse_functions.alt_depth("44,29,42", "0.477273,0,0.619048") == '21,0,26'


def test_gnomad_link():
    assert (
        parse_functions.gnomad_link("1:1581136:G:A")
        == '=HYPERLINK("http://gnomad.broadinstitute.org/variant/1-1581136-G-A","GNOMAD_link")'
    )


def test_ucsc_link():
    assert (
        parse_functions.ucsc_link("1:1581136:G:A")
        == '=HYPERLINK("http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&hgt.out3=10x&position=1:1581136","UCSC_link")'
    )


def test_add_c4r_exome_db():
    import pandas as pd

    test_record = pd.DataFrame.from_dict({"chr:pos:ref:alt": ["1:12854401:G:T"]})
    truth_record = pd.DataFrame.from_dict(
        {
            "chr:pos:ref:alt": ["1:12854401:G:T"],
            "Frequency_in_C4R": [3],
            "Seen_in_C4R_samples": ["207_120901A; 1469_CH0945; 1469_CH0950"],
        }
    )
    assert parse_functions.add_c4r_exome_db(test_record).equals(truth_record)


cano_refseq = 'MORN1/missense_variant/ENST00000378529/ENSG00000116151//9/11/ENST00000378529.3:c.790G>A/ENSP00000367790.3:p.Val264Met/264/350/missense_variant/probably_damaging(0.973)/deleterious(0.01)/PANTHER:PTHR23084&PANTHER:PTHR23084:SF123;MORN1/missense_variant/ENST00000378531/ENSG00000116151/YES/9/14/ENST00000378531.3:c.790G>A/ENSP00000367792.3:p.Val264Met/264/497/missense_variant/probably_damaging(0.939)/deleterious(0.02)/PANTHER:PTHR23084&PANTHER:PTHR23084:SF123;MORN1/non_coding_transcript_exon_variant/ENST00000606372/ENSG00000116151//8/12/ENST00000606372.1:n.872G>A///non_coding_transcript_exon_variant///;MORN1/non_coding_transcript_exon_variant/ENST00000607342/ENSG00000116151//4/4/ENST00000607342.1:n.533G>A///non_coding_transcript_exon_variant///;RP4-740C4.7/downstream_gene_variant/ENST00000607858/ENSG00000272420/YES/////downstream_gene_variant///;MORN1/missense_variant/NM_001301060.2/79906//9/11/NM_001301060.2:c.790G>A/NP_001287989.1:p.Val264Met/264/350/missense_variant/probably_damaging(0.973)/deleterious(0.01)/;MORN1/missense_variant/NM_024848.3/79906/YES/9/14/NM_024848.3:c.790G>A/NP_079124.1:p.Val264Met/264/497/missense_variant/probably_damaging(0.939)/deleterious(0.02)/'
refseq = 'ELMO1/5_prime_UTR_variant/ENST00000310758/ENSG00000155849/YES/1/22/ENST00000310758.4:c.-256_-253dup///5_prime_UTR_variant///;ELMO1/5_prime_UTR_variant/ENST00000445322/ENSG00000155849//1/5/ENST00000445322.1:c.-352_-349dup///5_prime_UTR_variant///;ELMO1/splice_region_variant&intron_variant/ENST00000448602/ENSG00000155849///ENST00000448602.1:c.-74+3_-74+6dup///splice_region_variant&intron_variant///;ELMO1/splice_region_variant&intron_variant/ENST00000453399/ENSG00000155849///ENST00000453399.1:c.-170+3_-170+6dup///splice_region_variant&intron_variant///;ELMO1/non_coding_transcript_exon_variant/ENST00000463390/ENSG00000155849//1/5/ENST00000463390.1:n.3_4insAAGT///non_coding_transcript_exon_variant///;ELMO1/splice_region_variant&intron_variant&non_coding_transcript_variant/ENST00000479447/ENSG00000155849///ENST00000479447.1:n.91+3_91+6dup///splice_region_variant&intron_variant&non_coding_transcript_variant///;ELMO1/splice_region_variant&intron_variant/NM_001206480.2/9844///NM_001206480.2:c.-74+3_-74+6dup///splice_region_variant&intron_variant///;ELMO1/5_prime_UTR_variant/NM_014800.11/9844//1/22/NM_014800.11:c.-256_-253dup///5_prime_UTR_variant///;/regulatory_region_variant/ENSR00000211027///////regulatory_region_variant///;/regulatory_region_variant/ENSR00001393429///////regulatory_region_variant///'
norefseq = 'ARAP2/splice_region_variant&intron_variant&non_coding_transcript_variant/ENST00000503225/ENSG00000047365///ENST00000503225.1:n.1470-8del///splice_region_variant&intron_variant&non_coding_transcript_variant///'

def test_parse_inf():
    # result = [refseq_HGVC, HGVCp, exon, protein, canonical, protein_domain, polyphen, sift]
    assert(
        parse_functions.parse_inf(cano_refseq) == ['NM_024848.3:c.790G>A', 'NP_079124.1:p.Val264Met', '9_14', '264_497', 'Yes', '', 'probably_damaging(0.939)', 'deleterious(0.02)']
    )
    assert(
        parse_functions.parse_inf(refseq) == ['NM_001206480.2:c.-74+3_-74+6dup', '', '', 'NA', 'No', '', '', '']
    )
    assert(
        parse_functions.parse_inf(norefseq) == ['NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA']
    )