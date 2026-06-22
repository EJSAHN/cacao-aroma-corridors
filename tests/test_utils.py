from cacao_aroma_pipeline.utils import coerce_integer, normalize_chromosome, slugify


def test_slugify():
    assert slugify("Cocoa Nacional aroma") == "cocoa_nacional_aroma"


def test_normalize_chromosome():
    assert normalize_chromosome("NC_030850.1_chromosome_1") == "1"
    assert normalize_chromosome("NC_030850.1") == "NC_030850.1"
    assert normalize_chromosome("chr10") == "10"


def test_coerce_integer():
    assert coerce_integer(" 1,234 ") == 1234
    assert coerce_integer("nan") is None
