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
