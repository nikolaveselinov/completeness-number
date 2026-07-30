from scripts.audit_census import audit


def test_public_census_and_degree_nine_certificates() -> None:
    result = audit()
    assert result["characteristics"] == 1352
    assert result["families"] == 31
    assert result["binary_characteristics"] == 127
    assert result["degree9_distribution"] == {14: 2, 16: 54}
    assert result["degree14_pairs"] == 14706
    assert result["degree13_zero_levels_each"] == 92
    assert result["degree9_manifest_files"] == 9
    assert result["detailed_rows"] == 1266
    assert result["degree8_rows"] == 30
    assert result["degree9_rows"] == 56


if __name__ == "__main__":
    test_public_census_and_degree_nine_certificates()
