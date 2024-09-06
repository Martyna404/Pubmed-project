from Bio import Entrez
import pandas as pd
import re
import time

Entrez.email = "martynapradela@gmail.com"

def clean_links(links):
    return [link.rstrip(').],') for link in links]

def get_pmids():
    total_pmids = []
    size = 1000  
    max_records = 9999  
    delay = 1  

    for start in range(0, max_records, size):
        stream = Entrez.esearch(db="pubmed", term="coronary heart disease risk factor", retstart=start, retmax=size, sort="relevance")
        results = Entrez.read(stream)
        pmids = results.get("IdList", [])
        total_pmids.extend(pmids)

        time.sleep(delay)

    return total_pmids

def get_data():
    pmids = get_pmids()
    stream = Entrez.efetch(db="pubmed", id=",".join(pmids), retmode="xml")
    results = Entrez.read(stream)
    stream.close()

    all_articles_data = []

    for article in results['PubmedArticle']:
        pmid = article['MedlineCitation']['PMID']

        e_location_id = article['PubmedData']['ArticleIdList']
        doi_url = 'DOI not found'
        for eloc in e_location_id:
            if eloc.attributes.get('IdType') == 'doi':
                doi_url = f"https://doi.org/{eloc}"

        title = article['MedlineCitation']['Article']['ArticleTitle']
        publication_type = ", ".join([ptype for ptype in article['MedlineCitation']['Article']['PublicationTypeList']])

        authors = ", ".join([
            f"{author.get('LastName', '')} {author.get('ForeName', '')}".strip() 
            for author in article['MedlineCitation']['Article'].get('AuthorList', [])
        ])
        
        mesh_terms = []
        for mesh_heading in article['MedlineCitation'].get('MeshHeadingList', []):
            descriptor = mesh_heading['DescriptorName']
            qualifiers = [qualifier for qualifier in mesh_heading.get('QualifierName', [])]
            mesh_term = descriptor
            if qualifiers:
                mesh_term += ": " + ", ".join(qualifiers)
            mesh_terms.append(mesh_term)

        mesh_terms_str = "; ".join(mesh_terms)
        country = article['MedlineCitation']['MedlineJournalInfo'].get('Country', 'Country not found')
        pub_date = article['MedlineCitation']['Article'].get('Journal', {}).get('JournalIssue', {}).get('PubDate', {})
        pub_date_str = f"{pub_date.get('Year', 'No Year')}"

        abstract_text = ''
        links_in_abstract = 'No Abstract'
        git_hub_links = 'No github links'  # Zainicjalizuj zmienną na początku pętli

        if 'Abstract' in article['MedlineCitation']['Article']:
            abstract = article['MedlineCitation']['Article']['Abstract']['AbstractText']
            if isinstance(abstract, list):
                abstract_text = " ".join(abstract)

            links = re.findall(r'https://\S+', abstract_text)
            cleaned_links = clean_links(links)
            links_in_abstract = ", ".join(cleaned_links) if cleaned_links else 'No links found'

            git_hub = re.findall(r'https://github\S+', abstract_text)
            cleaned_git_hub_links = clean_links(git_hub)
            git_hub_links = ", ".join(cleaned_git_hub_links) if cleaned_git_hub_links else 'No github links'

        article_data = {
            'PMID': pmid,
            'Abstract Links': links_in_abstract,
            'Github links': git_hub_links,
            'DOI URL': doi_url,
            'Title':  title,
            'Publication Year': publication_type,
            'Authors': authors,
            'Country': country,
            'Has Abstract': 'Yes' if abstract_text else 'No',
            'Publication Date': pub_date_str,
            'MeSH Terms': mesh_terms_str
        }

        all_articles_data.append(article_data)

    df = pd.DataFrame(all_articles_data)
    
    
    df.to_csv('coronaryheartdiseaseriskfactor.csv', index=False)
    
    print(df)
    return df


df = get_data()





""" 
tylko dla id=1

efetch
api1= https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=1&retmode=xml

esumarry
api2=https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=1

esearch
api3=https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=cancer&retmax=10&retmode=xml

"""
"""For PubMed, the valid sort orders are:
First Author
Journal
Last Author
Pub Date
Recently Added
Relevance
Title"""
