import os
import pytest

import pubmed_parser as pp


def test_parse_medline_xml():
    """
    Test parsing MEDLINE XML
    """
    expected_title = "Monitoring of bacteriological contamination and as"
    expected_abstract = "Two hundred and sixty nine beef, 230 sheep and 165"

    parsed_medline = pp.parse_medline_xml(os.path.join("data", "pubmed20n0014.xml.gz"))
    assert isinstance(parsed_medline, list)
    assert len(parsed_medline) == 30000, "Expect to have 30000 records"
    assert (
        len([p for p in parsed_medline if len(p["title"]) > 0]) == 30000
    ), "Expect every records to have title"
    assert parsed_medline[0]["title"][0:50] == expected_title
    assert parsed_medline[0]["issue"] == "50(2)"
    assert parsed_medline[0]["pages"] == "123-33"
    assert parsed_medline[0]["abstract"][0:50] == expected_abstract
    assert parsed_medline[0]["pmid"] == "399296"
    assert parsed_medline[0]["languages"] == "eng"
    assert parsed_medline[0]["vernacular_title"] == ""

    # Checking mesh terms on article with more possibilities
    assert parsed_medline[2]["mesh_terms"] == "D000818:Animals; D000906:Antibodies; D004283:Dog Diseases; " \
                                              "D004285:Dogs; D056890:Eukaryota; D005455:Fluorescent Antibody " \
                                              "Technique; D011528:Protozoan Infections; D011529:Protozoan " \
                                              "Infections, Animal"
    assert parsed_medline[2]["mesh_subheadings"] == "Q000032:analysis; Q000276:immunology"
    assert parsed_medline[2]["mesh_major_topics"] == "Q000032:analysis; Q000276:immunology; D005455:Fluorescent " \
                                                     "Antibody Technique; D011529:Protozoan Infections, Animal"
    assert parsed_medline[2]["mesh_full_terms"] == "D000818:Animals; D000906/Q000032:Antibodies/analysis; " \
                                                   "D004283/Q000276:Dog Diseases/immunology; D004285:Dogs; " \
                                                   "D056890/Q000276:Eukaryota/immunology; D005455:Fluorescent " \
                                                   "Antibody Technique; D011528/Q000276:Protozoan " \
                                                   "Infections/immunology; D011529:Protozoan Infections, Animal"


def test_parse_medline_grant_id():
    """
    Test parsing grants from MEDLINE XML
    """
    grants = pp.parse_medline_grant_id(os.path.join("data", "pubmed20n0014.xml.gz"))
    assert isinstance(grants, list)
    assert isinstance(grants[0], dict)
    assert grants[0]["pmid"] == "399300"
    assert grants[0]["grant_id"] == "HL17731"
    assert len(grants) == 484, "Expect number of grants in a given file to be 484"
