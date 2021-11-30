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
    #DP="84,72,69"
	#AB="0.238095,0,0.246377"
    assert parse_functions.alt_depth("84,72,69", "0.238095,0,0.246377") == [19.99998, 0.0, 17.000013]
    assert parse_functions.alt_depth("44,29,42", "0.477273,0,0.619048") == [21.000012, 0.0, 26.000016000000002]

def test_gnomad_link():
    assert parse_functions.gnomad_link("1:1581136:G:A") == "=HYPERLINK(\"http://gnomad.broadinstitute.org/variant/1-1581136-G-A\",\"GNOMAD_link\")"

def test_ucsc_link():
    assert parse_functions.ucsc_link("1:1581136:G:A") == '=HYPERLINK("http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&hgt.out3=10x&position=1:1581136","UCSC_link")'

def test_replacedelim():
    assert parse_functions.replace_comma("1,1,0") == '1/1/0'
    assert parse_functions.replace_comma("1,1,2") == '1/1/2'
    assert parse_functions.replace_comma("0.298,0.265,0.298") == '0.298/0.265/0.298'
    assert parse_functions.replace_comma("0.298,0.265,NA") == '0.298/0.265/NA'

def test_add_c4r_exome_db():
    import pandas as pd
    test_record = pd.DataFrame.from_dict({"chr:pos:ref:alt": ["1:12854401:G:T"]})
    truth_record = pd.DataFrame.from_dict({"chr:pos:ref:alt": ["1:12854401:G:T"], "Frequency_in_C4R":[3], "Seen_in_C4R_samples": ["207_120901A; 1469_CH0945; 1469_CH0950"]})
    assert parse_functions.add_c4r_exome_db(test_record).equals(truth_record)
    
#def test_parse_consequences():
#    consequence = 'CDK11B/missense_variant&splice_region_variant/ENST00000317673/ENSG00000248333//5/21/ENST00000317673.7:c.404C>T/ENSP00000463200.1:p.Thr135Met/135/783/missense_variant&splice_region_variant/benign(0.328)/deleterious_low_confidence(0)/;CDK11B/missense_variant&splice_region_variant/ENST00000340677/ENSG00000248333//4/20/ENST00000340677.5:c.355C>T/ENSP00000464016.1:p.Arg119Cys/119/772/missense_variant&splice_region_variant/possibly_damaging(0.703)/deleterious_low_confidence(0)/;CDK11B/missense_variant&splice_region_variant/ENST00000341832/ENSG00000248333//4/20/ENST00000341832.6:c.253C>T/ENSP00000463048.1:p.Arg85Cys/85/738/missense_variant&splice_region_variant/possibly_damaging(0.703)/deleterious_low_confidence(0)/;CDK11B/missense_variant&splice_region_variant/ENST00000407249/ENSG00000248333/YES/4/21/ENST00000407249.3:c.388C>T/ENSP00000464036.1:p.Arg130Cys/130/785/missense_variant&splice_region_variant/benign(0.415)/deleterious_low_confidence(0)/;CDK11B/upstream_gene_variant/ENST00000513088/ENSG00000248333//////upstream_gene_variant///;CDK11B/intron_variant/NM_001291345.2/984///NM_001291345.2:c.226-511C>T///intron_variant///;CDK11B/intron_variant/NM_001787.3/984/YES//NM_001787.3:c.232-511C>T///intron_variant///;CDK11B/intron_variant/NM_033486.3/984///NM_033486.3:c.232-511C>T///intron_variant///;CDK11B/intron_variant/NM_033487.3/984///NM_033487.3:c.-276-3771C>T///intron_variant///;CDK11B/intron_variant/NM_033489.3/984///NM_033489.3:c.130-511C>T///intron_variant///;CDK11B/intron_variant/NM_033490.3/984///NM_033490.3:c.-282-3265C>T///intron_variant///;/regulatory_region_variant/ENSR00000918659///////regulatory_region_variant///'

