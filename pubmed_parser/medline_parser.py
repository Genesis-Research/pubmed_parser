"""
Parsers for MEDLINE XML
"""
import re
import numpy as np
from itertools import chain
from collections import defaultdict
from pubmed_parser.utils import read_xml, stringify_children, month_or_day_formater

__all__ = ["parse_medline_xml", "parse_medline_grant_id"]


def parse_pmid(pubmed_article):
    """
    A function to parse PMID from a given Pubmed Article tree

    Parameters
    ----------
    pubmed_article: Element
        The lxml node pointing to a medline document

    Returns
    -------
    pmid: str
        A string of PubMed ID parsed from a given
    """
    medline = pubmed_article.find("MedlineCitation")
    if medline.find("PMID") is not None:
        pmid = medline.find("PMID").text
        return pmid
    else:
        article_ids = pubmed_article.find("PubmedData/ArticleIdList")
        if article_ids is not None:
            pmid = article_ids.find('ArticleId[@IdType="pmid"]')
            if pmid is not None:
                if pmid.text is not None:
                    pmid = pmid.text.strip()
                else:
                    pmid = ""
            else:
                pmid = ""
        else:
            pmid = ""
    return pmid


def parse_doi(pubmed_article):
    """
    A function to parse DOI from a given Pubmed Article tree

    Parameters
    ----------
    pubmed_article: Element
        The lxml node pointing to a medline document

    Returns
    -------
    doi: str
        A string of DOI parsed from a given ``pubmed_article``
    """
    medline = pubmed_article.find("MedlineCitation")
    article = medline.find("Article")
    elocation_ids = article.findall("ELocationID")

    if len(elocation_ids) > 0:
        for e in elocation_ids:
            doi = e.text.strip() or "" if e.attrib.get("EIdType", "") == "doi" else ""
    else:
        article_ids = pubmed_article.find("PubmedData/ArticleIdList")
        if article_ids is not None:
            doi = article_ids.find('ArticleId[@IdType="doi"]')
            doi = (
                (doi.text.strip() if doi.text is not None else "")
                if doi is not None
                else ""
            )
        else:
            doi = ""
    return doi


def parse_elocation_ids(pubmed_article):
    """
    A function to parse DOI from a given Pubmed Article tree

    Parameters
    ----------
    pubmed_article: Element
        The lxml node pointing to a medline document

    Returns
    -------
    doi: list|objects
        A list of ELocationID objects parsed from a given ``pubmed_article``. Added DOI!
    """
    medline = pubmed_article.find("MedlineCitation")
    article = medline.find("Article")
    elocation_ids_xml = article.findall("ELocationID")
    found_doi = False
    elocation_ids = list()

    if len(elocation_ids_xml) > 0:
        for e in elocation_ids_xml:
            elocation_id_type = e.attrib.get('EIdType')
            if elocation_id_type == "doi":
                found_doi = True
            elocation_ids.append({'type': elocation_id_type, 'ELocationID': e.text.strip()})

    if not found_doi:
        article_ids = pubmed_article.find("PubmedData/ArticleIdList")
        if article_ids is not None:
            doi = article_ids.find('ArticleId[@IdType="doi"]')
            if doi is not None:
                if doi.text.strip() != "":
                    elocation_ids.append({'type': 'doi', 'ELocationID': doi.text.strip()})

    return elocation_ids


def parse_mesh_terms(medline):
    """
    A function to parse MESH terms from article

    Parameters
    ----------
    medline: Element
        The lxml node pointing to a medline document

    Returns
    -------
    mesh_terms: str
        String of semi-colon ``;`` spearated MeSH (Medical Subject Headings)
        terms contained in the document.
    """
    if medline.find("MeshHeadingList") is not None:
        mesh = medline.find("MeshHeadingList")
        mesh_terms_list = [
            m.find("DescriptorName").attrib.get("UI", "")
            + ":"
            + m.find("DescriptorName").text
            for m in mesh.getchildren()
        ]
        mesh_terms = "; ".join(mesh_terms_list)
    else:
        mesh_terms = ""
    return mesh_terms


def parse_mesh_info(medline):
    """
    A function to parse MESH info from article, including MESH Headings, Subheadings, Major Topics, and Full Terms.

    Parameters
    ----------
    medline: Element
        The lxml node pointing to a medline document

    Returns
    -------
    dict_out: dict
        dictionary with keys including `mesh_terms`, `mesh_subheadings` `mesh_major_topics`, and `mesh_full_terms`

        mesh_terms: str
            String of semi-colon ``;`` spearated MeSH (Medical Subject Headings)
            terms contained in the document.
        mesh_subheadings: str
            String of semi-colon ``;`` separated MeSH (Medical Subject Headings)
            subheadings contained in the document.
        mesh_major_topics: str
            String of semi-colon ``;`` spearated MeSH (Medical Subject Headings)
            headings and subheadings that are designated as Major Topics contained in the document.
        mesh_full_terms: str
            String of semi-colon ``;`` separated MeSH (Medical Subject Headings)
            full terms, including MeSH Headings, Subheadings, and Major Topics (Depicted as * at end of text)
            contained in the document.
    """
    if medline.find("MeshHeadingList") is not None:
        mesh = medline.find("MeshHeadingList")
        mesh_terms_list = []
        mesh_subheadings_list = []
        mesh_major_topics_list = []
        mesh_full_terms_list = []
        for mesh_heading in mesh.getchildren():
            mesh_term = mesh_heading.find("DescriptorName")
            subheadings = mesh_heading.findall("QualifierName")

            mesh_terms_list.append(
                mesh_term.attrib.get("UI", "")
                + ":"
                + mesh_term.text
            )

            if len(subheadings) > 0:
                for subheading in subheadings:
                    subheading_text = subheading.attrib.get("UI", "") + ":" + subheading.text
                    if subheading_text not in mesh_subheadings_list:
                        mesh_subheadings_list.append(subheading_text)
                    if subheading.attrib.get("MajorTopicYN", "N") == "Y":
                        if subheading_text not in mesh_major_topics_list:
                            mesh_major_topics_list.append(subheading_text)

                    mesh_full_terms_list.append(
                        mesh_term.attrib.get("UI", "")
                        + "/"
                        + subheading.attrib.get("UI", "")
                        + ":"
                        + mesh_term.text
                        + "/"
                        + subheading.text
                        + ("*" if subheading.attrib.get("MajorTopicYN", "N") == "Y" else "")
                    )
            else:
                if mesh_term.attrib.get("MajorTopicYN", "N") == "Y":
                    mesh_major_topics_list.append(
                        mesh_term.attrib.get("UI", "")
                        + ":"
                        + mesh_term.text
                    )

                mesh_full_terms_list.append(
                    mesh_term.attrib.get("UI", "")
                    + ":"
                    + mesh_term.text
                    + ("*" if mesh_term.attrib.get("MajorTopicYN", "N") == "Y" else "")
                    )

        mesh_terms = "; ".join(mesh_terms_list)
        mesh_subheadings = "; ".join(mesh_subheadings_list)
        mesh_major_topics = "; ".join(mesh_major_topics_list)
        mesh_full_terms = "; ".join(mesh_full_terms_list)
    else:
        mesh_terms = ""
        mesh_subheadings = ""
        mesh_major_topics = ""
        mesh_full_terms = ""
    return {
        "mesh_terms": mesh_terms,
        "mesh_subheadings": mesh_subheadings,
        "mesh_major_topics": mesh_major_topics,
        "mesh_full_terms": mesh_full_terms
    }


def parse_publication_types(medline):
    """Parse Publication types from article

    Parameters
    ----------
    medline: Element
        The lxml node pointing to a medline document

    Returns
    -------
    publication_types: str
        String of semi-colon spearated publication types
    """
    publication_types = []
    publication_type_list = medline.find("Article/PublicationTypeList")
    if publication_type_list is not None:
        publication_type_list = publication_type_list.findall("PublicationType")
        for publication_type in publication_type_list:
            publication_types.append(
                publication_type.attrib.get("UI", "")
                + ":"
                + (publication_type.text.strip() or "")
            )
    publication_types = "; ".join(publication_types)
    return publication_types


def parse_keywords(medline):
    """Parse keywords from article, separated by ;

    Parameters
    ----------
    medline: Element
        The lxml node pointing to a medline document

    Returns
    -------
    keywords: str
        String of concatenated keywords.
    """
    keyword_list = medline.find("KeywordList")
    keywords = list()
    if keyword_list is not None:
        for k in keyword_list.findall("Keyword"):
            if k.text is not None:
                keywords.append(k.text)
        keywords = "; ".join(keywords)
    else:
        keywords = ""
    return keywords


def parse_chemical_list(medline):
    """Parse chemical list from article

    Parameters
    ----------
    medline: Element
        The lxml node pointing to a medline document

    Returns
    -------
    chemical_list: str
        String of semi-colon spearated chemical list
    """
    chemical_list = []
    chemicals = medline.find("ChemicalList")
    if chemicals is not None:
        for chemical in chemicals.findall("Chemical"):
            registry_number = chemical.find("RegistryNumber")
            substance_name = chemical.find("NameOfSubstance")
            chemical_list.append(
                (registry_number.text.strip() or "")
                + ":"
                + substance_name.attrib.get("UI", "")
                + ":"
                + (substance_name.text.strip() or "")
            )
    chemical_list = "; ".join(chemical_list)
    return chemical_list


def parse_other_id(medline):
    """Parse OtherID from article, each separated by ;

    Parameters
    ----------
    medline: Element
        The lxml node pointing to a medline document

    Returns
    -------
    other_id: str
        String of semi-colon separated Other IDs found in the document
    """
    pmc = ""
    other_id = list()
    oids = medline.findall("OtherID")
    if oids is not None:
        for oid in oids:
            if "PMC" in oid.text:
                pmc = oid.text
            else:
                other_id.append(oid.text)
        other_id = "; ".join(other_id)
    else:
        other_id = ""
    return {"pmc": pmc, "other_id": other_id}


def parse_journal_info(medline):
    """Parse MEDLINE journal information

    Parameters
    ----------
    medline: Element
        The lxml node pointing to a medline document

    Returns
    -------
    dict_out: dict
        dictionary with keys including `medline_ta`, `nlm_unique_id`,
        `issn_linking` and `country`
    """
    journal_info = medline.find("MedlineJournalInfo")
    if journal_info is not None:
        if journal_info.find("MedlineTA") is not None:
            medline_ta = (
                journal_info.find("MedlineTA").text or ""
            )  # equivalent to Journal name
        else:
            medline_ta = ""
        if journal_info.find("NlmUniqueID") is not None:
            nlm_unique_id = journal_info.find("NlmUniqueID").text or ""
        else:
            nlm_unique_id = ""
        if journal_info.find("ISSNLinking") is not None:
            issn_linking = journal_info.find("ISSNLinking").text
        else:
            issn_linking = ""
        if journal_info.find("Country") is not None:
            country = journal_info.find("Country").text or ""
        else:
            country = ""
    else:
        medline_ta = ""
        nlm_unique_id = ""
        issn_linking = ""
        country = ""
    dict_info = {
        "medline_ta": medline_ta.strip(),
        "nlm_unique_id": nlm_unique_id,
        "issn_linking": issn_linking,
        "country": country,
    }
    return dict_info


def parse_grant_id(pubmed_article):
    """Parse Grant ID and related information from a given MEDLINE tree

    Parameters
    ----------
    pubmed_article: Element
        The lxml node pointing to a medline document

    Returns
    -------
    grant_list: list
        List of grants acknowledged in the publications. Each
        entry in the dictionary contains the PubMed ID,
        grant ID, grant acronym, country, and agency.
    """
    medline = pubmed_article.find("MedlineCitation")
    article = medline.find("Article")
    pmid = parse_pmid(pubmed_article)

    grants = article.find("GrantList")
    grant_list = list()
    if grants is not None:
        grants_list = grants.getchildren()
        for grant in grants_list:
            grant_country = grant.find("Country")
            if grant_country is not None:
                country = grant_country.text
            else:
                country = ""
            grant_agency = grant.find("Agency")
            if grant_agency is not None:
                agency = grant_agency.text
            else:
                agency = ""
            grant_acronym = grant.find("Acronym")
            if grant_acronym is not None:
                acronym = grant_acronym.text
            else:
                acronym = ""
            grant_id = grant.find("GrantID")
            if grant_id is not None:
                gid = grant_id.text
            else:
                gid = ""
            grant_dict = {
                "pmid": pmid,
                "grant_id": gid,
                "grant_acronym": acronym,
                "country": country,
                "agency": agency,
            }
            grant_list.append(grant_dict)
    return grant_list


def parse_author_affiliation(medline):
    """Parse MEDLINE authors and their corresponding affiliations

    Parameters
    ----------
    medline: Element
        The lxml node pointing to a medline document

    Returns
    -------
    authors: list
        List of authors and their corresponding affiliation in dictionary format
    """
    authors = []
    article = medline.find("Article")
    if article is not None:
        # TODO: Make this work for books with multiple lists.
        if article.find("AuthorList") is not None:
            for author_list in article.findall("AuthorList"):
                authors_type = author_list.attrib.get("Type", "authors")
                authors_list = author_list.findall("Author")
                for author in authors_list:
                    if author.find("ForeName") is not None:
                        forename = (author.find("ForeName").text or "").strip() or ""
                    else:
                        forename = ""
                    if author.find("Initials") is not None:
                        initials = (author.find("Initials").text or "").strip() or ""
                    else:
                        initials = ""
                    if author.find("LastName") is not None:
                        lastname = (author.find("LastName").text or "").strip() or ""
                    else:
                        lastname = ""
                    if author.find("Identifier") is not None:
                        identifier_xml = author.find("Identifier")
                        identifier_type = (identifier_xml.attrib.get("Source", "") or "").strip() or ""
                        identifier = (identifier_xml.text or "").strip() or ""
                    else:
                        identifier_type = ""
                        identifier = ""
                    if author.find("Suffix") is not None:
                        suffix = (author.find("Suffix").text or "").strip() or ""
                    else:
                        suffix = ""
                    if author.find("CollectiveName") is not None:
                        corporate = (author.find("CollectiveName").text or "").strip() or ""
                    else:
                        corporate = ""
                    if author.find("AffiliationInfo/Affiliation") is not None:
                        affiliation = author.find("AffiliationInfo/Affiliation").text or ""
                        affiliation = affiliation.replace(
                            "For a full list of the authors' affiliations please see the Acknowledgements section.",
                            "",
                        )
                    else:
                        affiliation = ""
                    authors.append(
                        {
                            "lastname": lastname,
                            "forename": forename,
                            "initials": initials,
                            "suffix": suffix,
                            "author_type": authors_type,
                            "identifier_type": identifier_type,
                            "identifier": identifier,
                            "corporate": corporate,
                            "affiliation": affiliation
                        }
                    )
    return authors


def get_date_info(date_xml, year_info_only, parse_time=False):
    """Extract Date information from an Article in the Medline dataset.

    Parameters
    ----------
    date_xml: Element
        Any Date field in the Medline dataset
    year_info_only: bool
        if True, this tool will only attempt to extract year information from PubDate.
        if False, an attempt will be made to harvest all available PubDate information.
        If only year and month information is available, this will yield a date of
        the form 'YYYY-MM'. If year, month and day information is available,
        a date of the form 'YYYY-MM-DD' will be returned.
    parse_time: bool
        if True, this tool will also extract hour and minute info from the date.
        If only year and month information is available, this will yield a date of
        the form 'YYYY-MM'. If year, month and day information is available,
        a date of the form 'YYYY-MM-DD' will be returned. And so on.

    Returns
    -------
    Date: str
        Date extracted from an article.
        Note: If year_info_only is False and a month could not be
        extracted this falls back to year automatically.
    """

    month = None
    day = None
    hour = None
    minute = None
    if date_xml.find("Year") is not None:
        year = date_xml.find("Year").text
        if not year_info_only:
            if date_xml.find("Month") is not None:
                month = month_or_day_formater(date_xml.find("Month").text)
                if date_xml.find("Day") is not None:
                    day = month_or_day_formater(date_xml.find("Day").text)
                    if date_xml.find("Hour") is not None:
                        hour = month_or_day_formater(date_xml.find("Hour").text)
                        if date_xml.find("Minute") is not None:
                            minute = month_or_day_formater(date_xml.find("Minute").text)
    elif date_xml.find("MedlineDate") is not None:
        year_text = date_xml.find("MedlineDate").text
        year = re.findall(r"\d{4}", year_text)
        if len(year) >= 1:
            year = year[0]
        else:
            year = ""
    else:
        year = ""

    if year_info_only or month is None:
        return year
    else:
        date = "-".join(str(x) for x in filter(None, [year, month, day]))
        if parse_time and hour is not None:
            time = ":".join(str(x) for x in filter(None, [hour, minute]))
            date = date + ' ' + time
        return date


def date_extractor(journal, year_info_only):
    """Extract PubDate information from an Article in the Medline dataset.

    Parameters
    ----------
    journal: Element
        The 'Journal' field in the Medline dataset
    year_info_only: bool
        if True, this tool will only attempt to extract year information from PubDate.
        if False, an attempt will be made to harvest all available PubDate information.
        If only year and month information is available, this will yield a date of
        the form 'YYYY-MM'. If year, month and day information is available,
        a date of the form 'YYYY-MM-DD' will be returned.

    Returns
    -------
    PubDate: str
        PubDate extracted from an article.
        Note: If year_info_only is False and a month could not be
        extracted this falls back to year automatically.
    """
    issue = journal.xpath("JournalIssue")[0]
    issue_date = issue.find("PubDate")

    return get_date_info(issue_date, year_info_only)


def parse_references(pubmed_article, reference_list):
    """Parse references from Pubmed Article

    Parameter
    ---------
    pubmed_article: Element
        The lxml element pointing to a medline document

    reference_list: bool
        if it is True, return a list of dictionary
        if it is False return a string of PMIDs seprated by semicolon ';'

    Return
    ------
    references: (list, str)
        if 'reference_list' is set to True, return a list of dictionary
        if 'reference_list' is set to False return a string of PMIDs seprated by semicolon ';'
    """
    references = []
    reference_list_data = pubmed_article.find("PubmedData/ReferenceList")
    if reference_list_data is not None:
        for ref in reference_list_data.findall("Reference"):
            citation = ref.find("Citation")
            if citation is not None:
                if citation.text is not None:
                    citation = citation.text.strip()
                else:
                    citation = ""
            else:
                citation = ""
            article_ids = ref.find("ArticleIdList")
            pmid = (
                article_ids.find('ArticleId[@IdType="pubmed"]')
                if article_ids is not None
                else None
            )
            if pmid is not None:
                if pmid.text is not None:
                    pmid = pmid.text.strip()
                else:
                    pmid = ""
            else:
                pmid = ""
            references.append({"citation": citation, "pmid": pmid})

    if reference_list:
        return references
    else:
        references = ";".join(
            [ref["pmid"] for ref in references if ref["pmid"] != ""]
        )
        return references


def parse_coi_statement(medline):
    """Parse COI Statement from article

    Parameters
    ----------
    medline: Element
        The lxml node pointing to a medline document

    Returns
    -------
    coi_statement: str
        Returns the COI Statement, or "" if doesn't exist.
    """
    coi_statement = ""
    if medline.find("CoiStatement") is not None:
        coi_statement = (medline.find("CoiStatement").text or "").strip() or ""
    return coi_statement


def parse_pubmed_history_dates(pubmed_article, year_info_only, parse_time):
    """Parse history dates from article

    Parameters
    ----------
    pubmed_article: Element
        The lxml element pointing to a medline document
    year_info_only: bool
        if True, this tool will only attempt to extract year information from PubDate.
        if False, an attempt will be made to harvest all available PubDate information.
        If only year and month information is available, this will yield a date of
        the form 'YYYY-MM'. If year, month and day information is available,
        a date of the form 'YYYY-MM-DD' will be returned.
    parse_time: bool
        if True, parse the hour and minute from the date and add to the string as "YYYY/MM/DD HH:MM"
        if False, just parse the date info, "YYYY/MM/DD".
        This requires year_info_only to be True

    Returns
    -------
    dates: (list, dict)
        Return a list of dictionary. Status and Date.
    """
    dates = list()
    if pubmed_article.find("PubmedData") is not None:
        if pubmed_article.find("PubmedData/History") is not None:
            pubmed_history = pubmed_article.find("PubmedData/History")
            for pubmed_pubdate in pubmed_history.findall("PubMedPubDate"):
                pub_status = pubmed_pubdate.attrib.get("PubStatus")
                pub_date = get_date_info(pubmed_pubdate, year_info_only, parse_time)
                dates.append({'status': pub_status, 'date': pub_date})
    return dates


def parse_completion_date(medline, year_info_only):
    """Parse Completion date from article

    Parameters
    ----------
    medline: Element
        The lxml node pointing to a medline document
    year_info_only: bool
        if True, this tool will only attempt to extract year information from PubDate.
        if False, an attempt will be made to harvest all available PubDate information.
        If only year and month information is available, this will yield a date of
        the form 'YYYY-MM'. If year, month and day information is available,
        a date of the form 'YYYY-MM-DD' will be returned.

    Returns
    -------
    date_completion: str
        String of date completed as "YYYY-MM-DD"
    """

    date_completion = ""

    if medline.find("DateCompleted") is not None:
        date_completion_xml = medline.find("DateCompleted")
        date_completion = get_date_info(date_completion_xml, year_info_only)

    return date_completion


def parse_modification_date(medline, year_info_only):
    """Parse Modification date from article

    Parameters
    ----------
    medline: Element
        The lxml node pointing to a medline document
    year_info_only: bool
        if True, this tool will only attempt to extract year information from PubDate.
        if False, an attempt will be made to harvest all available PubDate information.
        If only year and month information is available, this will yield a date of
        the form 'YYYY-MM'. If year, month and day information is available,
        a date of the form 'YYYY-MM-DD' will be returned.

    Returns
    -------
    date_modification: str
        String of last modified date as "YYYY-MM-DD"
    """

    date_modification = ""

    if medline.find("DateRevised") is not None:
        date_revised_xml = medline.find("DateRevised")
        date_modification = get_date_info(date_revised_xml, year_info_only)

    return date_modification


def parse_investigators_list(medline):
    """Parse MEDLINE investigators

    Parameters
    ----------
    medline: Element
        The lxml node pointing to a medline document

    Returns
    -------
    investigators: list
        List of investigators in dictionary format
    """
    investigators = []
    investigator_list = medline.find("InvestigatorList")
    if investigator_list is not None:
        investigators_list = investigator_list.findall("Investigator")
        for investigator in investigators_list:
            if investigator.find("ForeName") is not None:
                forename = (investigator.find("ForeName").text or "").strip() or ""
            else:
                forename = ""
            if investigator.find("Initials") is not None:
                initials = (investigator.find("Initials").text or "").strip() or ""
            else:
                initials = ""
            if investigator.find("LastName") is not None:
                lastname = (investigator.find("LastName").text or "").strip() or ""
            else:
                lastname = ""
            if investigator.find("Suffix") is not None:
                suffix = (investigator.find("Suffix").text or "").strip() or ""
            else:
                suffix = ""
            investigators.append(
                {
                    "lastname": lastname,
                    "forename": forename,
                    "initials": initials,
                    "suffix": suffix
                }
            )
    return investigators


def parse_databank_list(medline):
    """
    A function to parse databank list from a given Pubmed Article tree

    Parameters
    ----------
    medline: Element
        The lxml node pointing to a medline document

    Returns
    -------
    databank_list: list|objects
        A list of Databank objects with name and accession numbers parsed from a given ``pubmed_article``.
    """
    article = medline.find("Article")
    databank_list_xml = article.find("DataBankList")
    databank_list = list()
    if databank_list_xml:
        for databank_xml in databank_list_xml.findall("DataBank"):
            databank_name = databank_xml.find("DataBankName")
            accession_number_list_xml = databank_xml.find("AccessionNumberList")
            accession_number_list = list()
            if accession_number_list_xml:
                accession_number_list = [
                    accession_number.text
                    for accession_number in accession_number_list_xml.findall("AccessionNumber")
                ]
            databank_list.append({"name": databank_name, "accession_numbers": accession_number_list})

    return databank_list


def parse_personal_subject_names_list(medline):
    """Parse personal subject names from medline

    Parameters
    ----------
    medline: Element
        The lxml node pointing to a medline document

    Returns
    -------
    personal_subject_names: list
        List of personal subject names in dictionary format
    """
    personal_subject_names = []
    personal_subject_names_xml = medline.find("PersonalNameSubjectList")
    if personal_subject_names_xml is not None:
        for personal_subject_name in personal_subject_names_xml.findall("PersonalNameSubject"):
            if personal_subject_name.find("ForeName") is not None:
                forename = (personal_subject_name.find("ForeName").text or "").strip() or ""
            else:
                forename = ""
            if personal_subject_name.find("Initials") is not None:
                initials = (personal_subject_name.find("Initials").text or "").strip() or ""
            else:
                initials = ""
            if personal_subject_name.find("LastName") is not None:
                lastname = (personal_subject_name.find("LastName").text or "").strip() or ""
            else:
                lastname = ""
            if personal_subject_name.find("Suffix") is not None:
                suffix = (personal_subject_name.find("Suffix").text or "").strip() or ""
            else:
                suffix = ""
            personal_subject_name.append(
                {
                    "lastname": lastname,
                    "forename": forename,
                    "initials": initials,
                    "suffix": suffix
                }
            )
    return personal_subject_names


def parse_supplementary_concepts_list(medline):
    """Parse supplementary concepts from medline

    Parameters
    ----------
    medline: Element
        The lxml node pointing to a medline document

    Returns
    -------
    supplementary_concepts: list
        List of supplementary concepts in dictionary format
    """
    supplementary_concepts = []
    supplementary_concepts_xml = medline.find("SupplMeshList")
    if supplementary_concepts_xml is not None:
        supplementary_concepts = [
            {
                "type": supplementary_concept.attrib.get("Type"),
                "UI": supplementary_concept.attrib.get("UI"),
                "supplementary_concept": supplementary_concept.text
            }
            for supplementary_concept in supplementary_concepts_xml.findall("SupplMeshName")
        ]

    return supplementary_concepts


def parse_article_info(
    pubmed_article, year_info_only, nlm_category, author_list, reference_list, parse_time, history_dates_list,
        investigator_list, elocation_ids_list, databank_list, personal_subject_names_list, supplementary_concepts_list,
        grant_ids_list
):
    """Parse article nodes from Medline dataset

    Parameters
    ----------
    pubmed_article: Element
        The lxml element pointing to a medline document
    year_info_only: bool
        see more details in date_extractor()
    nlm_category: bool
        see more details in parse_medline_xml()
    author_list: bool
        if True, return output as list, else string concatenated by ';'
    reference_list: bool
        if True, parse reference list as an output
    parse_time: bool
        if True, parse time from history dates and add to history dates output.
    history_dates_list: bool
        if True, return output as list, else string concatenated by ';'
    investigator_list: bool
        if True, return output as list, else string concatenated by ';'
    elocation_ids_list: bool
        if True, return output as list, else string concatenated by ';'
    databank_list: bool
        if True, return output as list, else string concatenated by ';'
    personal_subject_names_list: bool
        if True, return output as list, else string concatenated by ';'
    supplementary_concepts_list: bool
        if True, return output as list, else string concatenated by ';'
    grant_ids_list: bool
        if True, return output as list, else string concatenated by ';'

    Returns
    -------
    article: dict
        Dictionary containing information about the article, including
        `title`, `abstract`, `journal`, `authors`, `affiliations`, `pubdate`,
        `pmid`, `other_id`, `mesh_terms`, `mesh_subheadings`, `mesh_major_topics`,
        `mesh_full_terms`, `pages`, `issue`, and `keywords`. The field
        `delete` is always `False` because this function parses
        articles that by definition are not deleted.
    """
    medline = pubmed_article.find("MedlineCitation")
    article = medline.find("Article")

    if article.find("ArticleTitle") is not None:
        title = stringify_children(article.find("ArticleTitle")).strip() or ""
    else:
        title = ""

    if article.find("Journal/JournalIssue/Volume") is not None:
        volume = article.find("Journal/JournalIssue/Volume").text or ""
    else:
        volume = ""

    if article.find("Language") is not None:
        languages = ";".join([language.text for language in article.findall("Language")])
    else:
        languages = ""

    if article.find("VernacularTitle") is not None:
        vernacular_title = stringify_children(article.find("VernacularTitle")).strip() or ""
    else:
        vernacular_title = ""

    if article.find("Journal/JournalIssue/Issue") is not None:
        issue = article.find("Journal/JournalIssue/Issue").text or ""
    else:
        issue = ""

    short_issue = issue

    if volume == "":
        issue = ""
    else:
        issue = f"{volume}({issue})"

    if article.find("Pagination/MedlinePgn") is not None:
        pages = article.find("Pagination/MedlinePgn").text or ""
    else:
        pages = ""

    category = "NlmCategory" if nlm_category else "Label"
    if article.find("Abstract/AbstractText") is not None:
        # parsing structured abstract
        if len(article.findall("Abstract/AbstractText")) > 1:
            abstract_list = list()
            for abstract in article.findall("Abstract/AbstractText"):
                section = abstract.attrib.get(category, "")
                if section != "UNASSIGNED":
                    abstract_list.append("\n")
                    abstract_list.append(abstract.attrib.get(category, ""))
                section_text = stringify_children(abstract).strip()
                abstract_list.append(section_text)
            abstract = "\n".join(abstract_list).strip()
        else:
            abstract = (
                stringify_children(article.find("Abstract/AbstractText")).strip() or ""
            )
    elif article.find("Abstract") is not None:
        abstract = stringify_children(article.find("Abstract")).strip() or ""
    else:
        abstract = ""

    authors_dict = parse_author_affiliation(medline)
    if not author_list:
        affiliations = ";".join(
            [
                author.get("affiliation", "")
                for author in authors_dict
                if author.get("affiliation", "") != ""
            ]
        )
        authors = ";".join(
            [
                author.get("lastname", "") + "|" + author.get("forename",   "") + "|" +
                author.get("initials",  "") + "|" + author.get("suffix", "") + "|" + author.get("author_type", "") +
                "|" + author.get("identifier_type", "") + "|" + author.get("identifier", "") + "|" +
                author.get("corporate", "")
                for author in authors_dict
            ]
        )
    else:
        authors = authors_dict

    investigators_dict = parse_investigators_list(medline)
    if not investigator_list:
        investigators = ";".join(
            [
                investigator.get("lastname", "") + "|" + investigator.get("forename", "") + "|" +
                investigator.get("initials", "") + "|" + investigator.get("suffix", "")
                for investigator in investigators_dict
            ]
        )
    else:
        investigators = investigators_dict

    history_dates = parse_pubmed_history_dates(pubmed_article, year_info_only, parse_time)
    if not history_dates_list:
        history_dates = ";".join(
            [
                history_date.get("status", "") + "|" + history_date.get("date", "")
                for history_date in history_dates
            ]
        )

    elocation_ids = parse_elocation_ids(pubmed_article)
    if not elocation_ids_list:
        elocation_ids = ";".join(
            [
                elocation_id.get("type", "") + "|" + elocation_id.get("ELocationID", "")
                for elocation_id in elocation_ids
            ]
        )

    databanks = parse_databank_list(medline)
    if not databank_list:
        databanks_info = list()
        for databank in databanks:
            databank_name = databank.get("name")
            accession_numbers = databank.get("accession_numbers")
            if len(accession_numbers) > 0:
                for accession_number in databank.get("accession_numbers"):
                    databank_term = databank_name + "/" + accession_number
                    databanks_info.append(databank_term)
            else:
                databanks_info.append(databank_name)
        databanks_info = ";".join(databanks_info)
    else:
        databanks_info = databanks

    personal_subject_names = parse_personal_subject_names_list(medline)
    if not personal_subject_names_list:
        personal_subject_names = ";".join(
            [
                personal_subject_name.get("lastname", "") + "|" + personal_subject_name.get("forename", "") + "|" +
                personal_subject_name.get("initials", "") + "|" + personal_subject_name.get("suffix", "")
                for personal_subject_name in personal_subject_names
            ]
        )

    supplementary_concepts = parse_supplementary_concepts_list(medline)
    if not supplementary_concepts_list:
        supplementary_concepts = ";".join(
            [
                supplementary_concept.get("type") + "|" + supplementary_concept.get("UI")
                + "|" + supplementary_concept.get("supplementary_concept")
                for supplementary_concept in supplementary_concepts
            ]
        )

    publication_status = pubmed_article.find("PubmedData/PublicationStatus").text

    grant_ids = parse_grant_id(pubmed_article)
    if not grant_ids_list:
        grant_ids = ";".join(
            [
                grant_id.get("grant_id") + "|" + grant_id.get("grant_acronym")
                + "|" + grant_id.get("country") + "|" + grant_id.get("agency")
                for grant_id in grant_ids
            ]
        )

    journal = article.find("Journal")
    journal_name = " ".join(journal.xpath("Title/text()"))

    pmid = parse_pmid(pubmed_article)
    doi = parse_doi(pubmed_article)
    references = parse_references(pubmed_article, reference_list)
    pubdate = date_extractor(journal, year_info_only)
    mesh_results = parse_mesh_info(medline)
    publication_types = parse_publication_types(medline)
    chemical_list = parse_chemical_list(medline)
    keywords = parse_keywords(medline)
    other_id_dict = parse_other_id(medline)
    journal_info_dict = parse_journal_info(medline)
    coi_statement = parse_coi_statement(medline)
    completion_date = parse_completion_date(medline, year_info_only)
    modification_date = parse_modification_date(medline, year_info_only)

    dict_out = {
        "title": title,
        "issue": issue,
        "short_issue": short_issue,
        "volume": volume,
        "pages": pages,
        "abstract": abstract,
        "journal": journal_name,
        "authors": authors,
        "pubdate": pubdate,
        "pmid": pmid,
        "mesh_terms": mesh_results.get("mesh_terms"),
        "mesh_subheadings": mesh_results.get("mesh_subheadings"),
        "mesh_major_topics": mesh_results.get("mesh_major_topics"),
        "mesh_full_terms": mesh_results.get("mesh_full_terms"),
        "publication_types": publication_types,
        "chemical_list": chemical_list,
        "keywords": keywords,
        "doi": doi,
        "references": references,
        "delete": False,
        "languages": languages,
        "vernacular_title": vernacular_title,
        "coi_statement": coi_statement,
        "completion_date": completion_date,
        "modification_date": modification_date,
        "history_dates": history_dates,
        "investigators": investigators,
        "elocation_ids": elocation_ids,
        "databanks": databanks_info,
        "personal_subject_names": personal_subject_names,
        "supplementary_concepts": supplementary_concepts,
        "publication_status": publication_status,
        "grant_ids": grant_ids
    }
    if not author_list:
        dict_out.update({"affiliations": affiliations})

    dict_out.update(other_id_dict)
    dict_out.update(journal_info_dict)
    return dict_out


def parse_medline_xml(
    path,
    year_info_only=True,
    nlm_category=False,
    author_list=False,
    reference_list=False,
    parse_time=False,
    history_dates_list=False,
    investigator_list=False,
    elocation_ids_list=False,
    databanks_list=False,
    personal_subject_names_list=False,
    supplementary_concepts_list=False,
    grant_ids_list=False
):
    """Parse XML file from Medline XML format available at
    ftp://ftp.nlm.nih.gov/nlmdata/.medleasebaseline/gz/

    Parameters
    ----------
    path: str
        The path
    year_info_only: bool
        if True, this tool will only attempt to extract year information from PubDate.
        if False, an attempt will be made to harvest all available PubDate information.
        If only year and month information is available, this will yield a date of
        the form 'YYYY-MM'. If year, month and day information is available,
        a date of the form 'YYYY-MM-DD' will be returned.
        NOTE: the resolution of PubDate information in the Medline(R) database varies
        between articles.
        default: True
    nlm_category: bool
        if True, this will parse structured abstract where each section if original Label
        if False, this will parse structured abstract where each section will be assigned to
        NLM category of each sections
        default: False
    author_list: bool
        if True, return parsed author output as a list of authors
        if False, return parsed author output as a string of authors concatenated with ``;``
        default: False
    reference_list: bool
        if True, parse reference list as an output
        if False, return string of PMIDs concatenated with ;
        default: False
    parse_time: bool
        if True, parse the time from history dates.
        if False, do not parse the time from the history dates.
        default: False
    history_dates_list: bool
        if True, return parsed history dates output as a list of dates
        if False, return parsed history output as a string of dates concatenated with ``;``
        default: False
    investigator_list: bool
        if True, return parsed investigator output as a list of investigators
        if False, return parsed investigator output as a string of investigators concatenated with ``;``
        default: False
    elocation_ids_list: bool
        if True, return parsed elocation_ids output as a list of elocation_id objects with type and elocation id.
        if False, return parsed elocation_ids output as a string of elocation_ids concatenated with ``;``
        default: False
    databanks_list: bool
        if True, return parsed databanks output as a list of databank objects with name and a list of accession numbers.
        if False, return parsed databanks output as a string of databanks concatenated with ``;``
        default: False
    personal_subject_names_list: bool
        if True, return parsed personal_subject_names output as a list of personal_subject_names objects.
        if False, return parsed personal_subject_names output as a string of subject_names concatenated with ``;``
        default: False
    supplementary_concepts_list: bool
        if True, return parsed supplementary_concepts output as a list of supplementary_concepts objects.
        if False, return parsed supplementary_concepts output as a string of concepts concatenated with ``;``
        default: False
    grant_ids_list: bool
        if True, return parsed grant_ids output as a list of grant_id objects.
        if False, return parsed grant_ids output as a string of grants concatenated with ``;`` (pmid arg not included)
        default: False

    Return
    ------
    article_list: list
        A list of dictionary containing information about articles in NLM format (see
        `parse_article_info`). Articles that have been deleted will be
        added with no information other than the field `delete` being `True`

    Examples
    --------
    >>> pubmed_parser.parse_medline_xml('data/pubmed20n0014.xml.gz')
    """
    tree = read_xml(path)
    medline_citations = tree.findall(".//MedlineCitationSet/MedlineCitation")
    if len(medline_citations) == 0:
        medline_citations = tree.findall(".//PubmedArticle")
    article_list = list(
        map(
            lambda m: parse_article_info(
                m, year_info_only, nlm_category, author_list, reference_list, parse_time, history_dates_list,
                investigator_list, elocation_ids_list, databanks_list, personal_subject_names_list,
                supplementary_concepts_list, grant_ids_list
            ),
            medline_citations,
        )
    )
    delete_citations = tree.findall(".//DeleteCitation/PMID")
    dict_delete = [
        {
            "title": np.nan,
            "abstract": np.nan,
            "journal": np.nan,
            "authors": np.nan,
            "affiliations": np.nan,
            "pubdate": np.nan,
            "pmid": p.text.strip(),
            "doi": np.nan,
            "other_id": np.nan,
            "pmc": np.nan,
            "mesh_terms": np.nan,
            "mesh_subheadings": np.nan,
            "mesh_major_topics": np.nan,
            "mesh_full_terms": np.nan,
            "keywords": np.nan,
            "publication_types": np.nan,
            "chemical_list": np.nan,
            "delete": True,
            "medline_ta": np.nan,
            "nlm_unique_id": np.nan,
            "issn_linking": np.nan,
            "country": np.nan,
            "references": np.nan,
            "issue": np.nan,
            "short_issue": np.nan,
            "volume": np.nan,
            "pages": np.nan,
            "languages": np.nan,
            "vernacular_title": np.nan,
            "coi_statement": np.nan,
            "completion_date": np.nan,
            "modification_date": np.nan,
            "history_dates": np.nan,
            "investigators": np.nan,
            "elocation_ids": np.nan,
            "databanks": np.nan,
            "personal_subject_names": np.nan,
            "supplementary_concepts": np.nan,
            "publication_status": np.nan,
            "grant_ids": np.nan
        }
        for p in delete_citations
    ]
    article_list.extend(dict_delete)
    return article_list


def parse_medline_grant_id(path):
    """Parse grant id from Medline XML file

    Parameters
    ----------
    path: str
        The path to the XML with the information

    Return
    ------
    grant_id_list: list
        A list of dictionaries contains the grants in a given path. Each dictionary
        has the keys of 'pmid', 'grant_id', 'grant_acronym', 'country', and 'agency'

    >>> pubmed_parser.parse_medline_grant_id('data/pubmed20n0014.xml.gz')
    [{
        'pmid': '399300',
        'grant_id': 'HL17731',
        'grant_acronym': 'HL',
        'country': 'United States',
        'agency': 'NHLBI NIH HHS'
    }, ...
    ]
    """
    tree = read_xml(path)
    medline_citations = tree.findall(".//MedlineCitationSet/MedlineCitation")
    if len(medline_citations) == 0:
        medline_citations = tree.findall(".//PubmedArticle")
    grant_id_list = list(map(parse_grant_id, medline_citations))
    grant_id_list = list(chain(*grant_id_list))  # flatten list
    return grant_id_list
