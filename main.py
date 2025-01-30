import io
import logging
import os
import re
import time
from urllib.error import HTTPError

import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import plotly.express as px
import streamlit as st
from Bio import Entrez
from matplotlib.colors import LinearSegmentedColormap
from wordcloud import WordCloud


Entrez.email = "martynapradela@gmail.com"

def setup_logger(term, start_year, end_year):
    log_folder = f"log-{term}"
    os.makedirs(log_folder, exist_ok=True)
    log_file = os.path.join(log_folder, "search.log")
    logger = logging.getLogger(term) 
    logger.setLevel(logging.INFO)  
    
    if logger.hasHandlers():
        logger.handlers.clear()

    file_handler = logging.FileHandler(log_file, encoding="utf-8")
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))

    
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    stream_handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)

    logger.info("### Rozpoczęcie wyszukiwania... ###")
    logger.info(f"Fraza wyszukiwania: {term}")
    logger.info(f"Zakres lat: {start_year} - {end_year}")
    
    
    return logger, time.time()  

def log_search_time(logger, start_time):
    total_time = time.time() - start_time
    logger.info(f"### Pobieranie zakończone ###")
    logger.info(f"Łączny czas pobierania danych: {total_time:.2f} sekundy.")


# Funkcja obsługująca zapytania z ponawianiem
def fetch_with_retry(fetch_function, *args, max_retries=5, logger=None, **kwargs):
    for attempt in range(1, max_retries + 1):
        try:
            return fetch_function(*args, **kwargs)
        except HTTPError as e:
            if e.code in [429, 500]:
                wait_time = 2 * attempt
                if logger:
                    logger.warning(f"Błąd HTTP {e.code}. Próba {attempt}/{max_retries}. Oczekiwanie {wait_time} sekund.")
                time.sleep(wait_time)
            else:
                if logger:
                    logger.error(f"Błąd HTTP {e.code}: {e.reason}.")
                raise
        except Exception as e:
            if logger:
                logger.error(f"Nieoczekiwany błąd: {e}. Próba {attempt}/{max_retries}.")
            time.sleep(2 * attempt)
        if attempt == max_retries:
            raise RuntimeError("Przekroczono maksymalną liczbę prób.")
  

def fetch_keywords(pmid):
    try:
        handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
        records = Entrez.read(handle)
        handle.close()

        article = records["PubmedArticle"][0]
        keywords_list = article["MedlineCitation"]["KeywordList"]

        if keywords_list:
            return ", ".join([kw for sublist in keywords_list for kw in sublist])
        return "Brak słów kluczowych"
    except Exception as e:
        return "Błąd pobierania"

def generate_wordcloud(text):
    wordcloud = WordCloud(width=800, height=400, background_color="white").generate(text)

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.imshow(wordcloud, interpolation="bilinear")
    ax.axis("off")
    return fig
        
def extract_github_links(abstract):
    if not abstract:
        return "No GitHub links"
    links = re.findall(r'https://github\.com/[^\s\)]+', abstract)
    cleaned_links = [link.rstrip("),.") for link in links]
    return ", ".join(cleaned_links) if cleaned_links else "No GitHub links"       



def fetch_articles_with_details(pmids, start_year, end_year, logger=None, progress_bar=None, total_steps=1, current_step=0):
    all_articles_data = []
    all_keywords = set()  
    github_repos_by_year = {}  
    
    batch_size = 1000  
    total_batches = len(pmids) // batch_size + (1 if len(pmids) % batch_size > 0 else 0)  

    logger.info(f"Rozpoczynanie pobierania szczegółowych danych dla {len(pmids)} PMIDs...")

    for i, start in enumerate(range(0, len(pmids), batch_size)):
        batch_pmids = pmids[start:start + batch_size]
        try:
            stream = Entrez.efetch(db="pubmed", id=",".join(batch_pmids), retmode="xml")
            results = Entrez.read(stream)
            stream.close()

            for article in results.get("PubmedArticle", []):
                try:
                    pmid = article["MedlineCitation"]["PMID"]
                    title = article["MedlineCitation"]["Article"]["ArticleTitle"]
                    country = article["MedlineCitation"]["MedlineJournalInfo"].get("Country", "Country not found")
                    pub_date = article["MedlineCitation"]["Article"].get("Journal", {}).get("JournalIssue", {}).get("PubDate", {})
                    pub_date_year = pub_date.get("Year", None)

                    if pub_date_year and (start_year <= int(pub_date_year) <= end_year):
                        authors = ", ".join([
                            f"{author.get('LastName', '')} {author.get('ForeName', '')}".strip()
                            for author in article["MedlineCitation"]["Article"].get("AuthorList", [])
                        ])

                        #Pobieranie słów kluczowych
                        keywords = set()
                        if "KeywordList" in article["MedlineCitation"]:
                            for kw_list in article["MedlineCitation"]["KeywordList"]:
                                if isinstance(kw_list, list):
                                    keywords.update(str(kw) for kw in kw_list)  
                                else:
                                    keywords.add(str(kw_list))  

                        #Pobieranie MeshHeadingList
                        if "MeshHeadingList" in article["MedlineCitation"]:
                            mesh_terms = article["MedlineCitation"]["MeshHeadingList"]
                            keywords.update(str(term["DescriptorName"]) for term in mesh_terms)

                        all_keywords.update(keywords)  
                        keywords_text = "; ".join(keywords) if keywords else "No Keywords"

                        
                        abstract_text = " ".join(article["MedlineCitation"]["Article"]["Abstract"]["AbstractText"]) if "Abstract" in article["MedlineCitation"]["Article"] else ""
                        github_links = extract_github_links(abstract_text)
                        has_github_links = "Yes" if "https://github.com" in github_links else "No"

                        
                        if has_github_links == "Yes":
                            github_repos_by_year[pub_date_year] = github_repos_by_year.get(pub_date_year, 0) + 1

                        article_data = {
                            "PMID": pmid,
                            "Title": title,
                            "Authors": authors,
                            "Country": country,
                            "Publication Date": pub_date_year,
                            "Abstract": abstract_text,
                            "Has GitHub Links": has_github_links,
                            "GitHub Links": github_links,
                            "Keywords": keywords_text,  
                        }
                        all_articles_data.append(article_data)

                except Exception as e:
                    logger.error(f"Błąd przetwarzania artykułu: {e}")
                    continue

            logger.info(f"Pobrano szczegóły dla {len(batch_pmids)} PMIDs.")

        except Exception as e:
            logger.error(f"Błąd podczas pobierania partii PMIDs: {e}")
            continue

        if progress_bar:
            progress = (current_step + (i + 1) / total_batches) / total_steps
            progress_bar.progress(progress)

    logger.info(f"Łącznie pobrano szczegóły dla {len(all_articles_data)} artykułów.")

   
    st.session_state["all_keywords"] = " ".join(all_keywords) if all_keywords else "No Keywords"
    st.session_state["github_repos_by_year"] = github_repos_by_year

    return all_articles_data



def get_pmids_with_full_pagination(term, start_year, end_year, logger=None, progress_bar=None, total_steps=1, current_step=0):
    total_pmids = []
    years = list(range(start_year, end_year + 1))
    total_years = len(years)

    logger.info("### Rozpoczynanie pobierania PMIDs ###")
    for i, year in enumerate(years):
        yearly_pmids = []
        try:
            year_query = f"{term} AND ({year}/01/01:{year}/12/31[dp])"

            def execute_query():
                return Entrez.read(Entrez.esearch(db="pubmed", term=year_query, retstart=0, retmax=9999))

            results = fetch_with_retry(execute_query, logger=logger)
            if results is None:
                logger.error(f"Nie udało się pobrać wyników dla roku {year}.")
                continue

            pmids = results.get("IdList", [])
            yearly_pmids.extend(pmids)

            logger.info(f"Pobrano {len(pmids)} wyników dla roku {year}.")

            # Jeśli liczba wyników przekracza 9999, dzielimy na miesiące
            if len(pmids) == 9999:
                logger.info(f"Dużo wyników dla roku {year}, dzielenie na miesiące...")
                monthly_pmids = get_pmids_for_months(term, year, logger)
                yearly_pmids.extend(monthly_pmids)

            total_pmids.extend(yearly_pmids)

        except Exception as e:
            logger.error(f"Błąd podczas pobierania wyników dla roku {year}: {e}")

        
        if progress_bar:
            progress = (current_step + (i + 1) / total_years) / total_steps
            progress_bar.progress(progress)

    unique_pmids = list(set(total_pmids))
    duplicate_count = len(total_pmids) - len(unique_pmids)

    logger.info("### Ogólne podsumowanie ###")
    logger.info(f"Łącznie pobrano {len(total_pmids)} PMIDs (z duplikatami).")
    logger.info(f"Liczba unikalnych PMIDs: {len(unique_pmids)}.")

    return unique_pmids, duplicate_count


def get_pmids_for_months(term, year, logger, delay=1):
    monthly_pmids = []
    for month in range(1, 13):
        try:
            month_query = f"{term} AND ({year}/{month:02d}/01:{year}/{month:02d}/31[dp])"

            def execute_query():
                return Entrez.read(Entrez.esearch(db="pubmed", term=month_query, retstart=0, retmax=9999))

            results = fetch_with_retry(execute_query, logger=logger)
            if results is None:
                logger.error(f"Nie udało się pobrać wyników dla miesiąca {year}-{month:02d}.")
                continue

            pmids = results.get("IdList", [])
            monthly_pmids.extend(pmids)
            logger.info(f"Pobrano {len(pmids)} wyników dla miesiąca {year}-{month:02d}.")

            # Jeśli liczba wyników przekracza 9999, dzielimy na dni
            if len(pmids) == 9999:
                logger.info(f"Dużo wyników dla miesiąca {year}-{month:02d}, dzielenie na dni...")
                daily_pmids = get_pmids_for_days(term, year, month, logger)
                monthly_pmids.extend(daily_pmids)

        except Exception as e:
            logger.error(f"Błąd podczas pobierania wyników dla miesiąca {year}-{month:02d}: {e}")
    return monthly_pmids


def get_pmids_for_days(term, year, month, logger, delay=1):
    daily_pmids = []
    for day in range(1, 32):  
        try:
            day_query = f"{term} AND ({year}/{month:02d}/{day:02d}:{year}/{month:02d}/{day:02d}[dp])"

            def execute_query():
                return Entrez.read(Entrez.esearch(db="pubmed", term=day_query, retstart=0, retmax=9999))

            results = fetch_with_retry(execute_query, logger=logger)
            if results is None:
                logger.error(f"Nie udało się pobrać wyników dla dnia {year}-{month:02d}-{day:02d}.")
                continue

            pmids = results.get("IdList", [])
            daily_pmids.extend(pmids)
            logger.info(f"Pobrano {len(pmids)} wyników dla dnia {year}-{month:02d}-{day:02d}.")
        except Exception as e:
            logger.error(f"Błąd podczas pobierania wyników dla dnia {year}-{month:02d}-{day:02d}: {e}")
    return daily_pmids

def create_interactive_bar_chart(data, output_file=None, logger=None):
    try:
        df = pd.DataFrame({
            "Year": list(data.keys()),
            "Number of Publications": list(data.values())
        })

        fig = px.bar(
            df,
            x="Year",
            y="Number of Publications",
            title="Liczba publikacji na przestrzeni lat",
            labels={"Year": "Rok", "Number of Publications": "Liczba publikacji"},
            hover_data={"Number of Publications": True},
        )

        
        fig.update_traces(marker=dict(line=dict(width=1, color="DarkSlateGrey")))
        fig.update_layout(
            hovermode="x",
            xaxis=dict(tickangle=45),
            template="plotly_white"
        )

        
        fig.show()
        if logger:
            logger.info("Wykres słupkowy został wygenerowany pomyślnie.")

        
        if output_file:
            if output_file.endswith(".html"):
                fig.write_html(output_file)
                if logger:
                    logger.info(f"Wykres zapisano jako plik HTML: {output_file}")
            elif output_file.endswith(".png"):
                fig.write_image(output_file)
                if logger:
                    logger.info(f"Wykres zapisano jako plik PNG: {output_file}")
            else:
                if logger:
                    logger.warning("Nieobsługiwany format pliku. Obsługiwane są: .html, .png")

    except Exception as e:
        if logger:
            logger.error(f"Błąd podczas tworzenia lub zapisu wykresu: {e}")
        raise

def create_map(dataframe):  
    country_counts = dataframe['Country'].value_counts().reset_index()
    country_counts.columns = ['Country', 'Count']
    world = gpd.read_file("C:/Users/Martyna/OneDrive/Pulpit/web scraping/ne_110m_admin_0_countries.shp")
    world = world.merge(country_counts, left_on='NAME', right_on='Country', how='left')
    world['Count'] = world['Count'].fillna(0)
    custom_cmap = LinearSegmentedColormap.from_list(
        "custom_cmap",
        ["grey", "yellow", "blue", "green"]
    )

    fig, ax = plt.subplots(1, 1, figsize=(15, 10))
    world.boundary.plot(ax=ax, linewidth=1)
    world.plot(
        column='Count',
        ax=ax,
        legend=True,
        cmap=custom_cmap,
        legend_kwds={'label': "Number of Publications", 'shrink': 0.7}
    )
    ax.set_title("Number of Publications by Country", fontsize=16)
    ax.set_axis_off()
    return fig 


##########STREAMLIT###############################

st.title("PubMed Data Viewer")
term = st.text_input("Wprowadź termin wyszukiwania:", value="heart")
start_year = st.number_input("Podaj początkowy rok:", value=2000, step=1)
end_year = st.number_input("Podaj końcowy rok:", value=2023, step=1)


if "articles_data" not in st.session_state:
    st.session_state["articles_data"] = None
if "publications_by_year" not in st.session_state:
    st.session_state["publications_by_year"] = None
if "github_repos_by_year" not in st.session_state:
    st.session_state["github_repos_by_year"] = None
if "all_keywords" not in st.session_state:
    st.session_state["all_keywords"] = ""

if st.button("Wyszukaj"):
    is_data_loaded = False
    with st.spinner("Pobieranie danych..."):
        try:
            logger, start_time = setup_logger(term, start_year, end_year)
            progress_bar = st.progress(0)

            
            pmids, duplicate_count = get_pmids_with_full_pagination(
                term, start_year, end_year, logger, progress_bar=progress_bar, total_steps=2, current_step=0
            )

            
            articles_data = fetch_articles_with_details(
                pmids, start_year, end_year, logger, progress_bar=progress_bar, total_steps=2, current_step=1
            )

            
            progress_bar.empty()

            
            log_search_time(logger, start_time)

            
            publications_by_year = {}
            github_repos_by_year = {}
            for article in articles_data:
                pub_date = article["Publication Date"]
                has_github = article["Has GitHub Links"] == "Yes"

                if pub_date.isdigit():
                    publications_by_year[pub_date] = publications_by_year.get(pub_date, 0) + 1
                    if has_github:
                        github_repos_by_year[pub_date] = github_repos_by_year.get(pub_date, 0) + 1

            
            st.session_state["articles_data"] = articles_data
            st.session_state["publications_by_year"] = publications_by_year
            st.session_state["github_repos_by_year"] = github_repos_by_year

            is_data_loaded = True

        except Exception as e:
            st.error(f"Wystąpił błąd: {e}")
            is_data_loaded = False


if st.session_state["articles_data"]:
    st.write("### Szczegółowe dane artykułów ###")
    df = pd.DataFrame(st.session_state["articles_data"])

    
    ROWS_PER_PAGE = 10000
    total_rows = len(df)
    total_pages = (total_rows // ROWS_PER_PAGE) + int(total_rows % ROWS_PER_PAGE > 0)
    page = st.number_input("Wybierz stronę:", min_value=1, max_value=total_pages, value=1, step=1, format="%d")

    start_row = (page - 1) * ROWS_PER_PAGE
    end_row = min(start_row + ROWS_PER_PAGE, total_rows)
    st.write(f"Wiersze {start_row + 1} do {end_row} z {total_rows}")
    st.dataframe(df.iloc[start_row:end_row])

    st.subheader(f"Wykres słupkowy dla: {term} w PubMed")
    publications_by_year = st.session_state["publications_by_year"]

    if publications_by_year:
        fig_bar = px.bar(
            x=list(publications_by_year.keys()),
            y=list(publications_by_year.values()),
            labels={"x": "Rok", "y": "Liczba publikacji"},
            title=f"Wykres słupkowy dla: {term} w PubMed"
        )

        
        fig_bar.update_xaxes(type="category")

        
        fig_bar.update_xaxes(tickmode="linear", dtick=1)

        
        fig_bar.update_traces(hovertemplate='Rok: %{x}<br>Liczba publikacji: %{y}<extra></extra>')

        st.plotly_chart(fig_bar, use_container_width=True)

        
        buf_bar = io.BytesIO()
        fig_bar.write_image(buf_bar, format="png", scale=2)
        buf_bar.seek(0)
        st.download_button("Pobierz wykres słupkowy jako PNG", buf_bar, "wykres_slupkowy.png", "image/png", key="download_bar_chart")

        
        st.subheader(f"Liczba znalezionych repozytoriów GitHub / rok")
        github_repos_by_year = st.session_state["github_repos_by_year"]

        if github_repos_by_year:
            fig_github = px.bar(
                x=list(map(int, github_repos_by_year.keys())),
                y=list(github_repos_by_year.values()),
                labels={"x": "Rok", "y": "Liczba repozytoriów GitHub"},
                title=f"Liczba znalezionych repozytoriów GitHub / rok"
            )
            fig_github.update_traces(hovertemplate='Rok: %{x}<br>Liczba repozytoriów: %{y}<extra></extra>')

            
            fig_github.update_layout(
                xaxis=dict(
                    tickmode="array",
                    tickvals=list(map(int, github_repos_by_year.keys())),
                    tickformat="d"
                )
            )

    st.plotly_chart(fig_github, use_container_width=True)

    st.write("### Mapa publikacji według krajów ###")
    map_figure = create_map(df)
    st.pyplot(map_figure)

    buf_map = io.BytesIO()
    map_figure.savefig(buf_map, format="png", dpi=300, bbox_inches="tight")
    buf_map.seek(0)
    st.download_button("Pobierz mapę jako PNG", buf_map, "mapa_publikacji.png", "image/png", key="download_map_chart")

       
    st.write("### TOP 5 krajów z największą liczbą publikacji")
    if "articles_data" in st.session_state and st.session_state["articles_data"]:
        df_countries = pd.DataFrame(st.session_state["articles_data"])
        if "Country" in df_countries.columns:
            country_counts = df_countries["Country"].value_counts().reset_index()
            country_counts.columns = ["Kraj", "Liczba publikacji"]
            top_5_countries = country_counts.head(5)

            table_html = top_5_countries.to_html(index=False, classes="styled-table")
            st.markdown(
                """
                <style>
                .styled-table {
                    width: 50%;
                    margin-left: auto;
                    margin-right: auto;
                    border-collapse: collapse;
                    text-align: center;
                    font-size: 16px;
                }
                .styled-table th, .styled-table td {
                    border: 1px solid #ddd;
                    padding: 10px;
                    text-align: center;
                }
                .styled-table th {
                    background-color: #f4f4f4;
                    font-weight: bold;
                }
                </style>
                """,
                unsafe_allow_html=True
            )

            st.markdown(table_html, unsafe_allow_html=True)

        else:
            st.write("Brak danych o krajach.")

    else:
        st.write("Brak dostępnych danych do analizy.")


    
    if "all_keywords" in st.session_state and st.session_state["all_keywords"] != "No Keywords":
        st.write("### Chmura słów kluczowych ###")

        keywords_text = st.session_state["all_keywords"]  

        wordcloud_fig = generate_wordcloud(keywords_text)
        st.pyplot(wordcloud_fig)

        
        buf_wordcloud = io.BytesIO()
        wordcloud_fig.savefig(buf_wordcloud, format="png", dpi=300, bbox_inches="tight")
        buf_wordcloud.seek(0)
        st.download_button("Pobierz chmurę słów jako PNG", buf_wordcloud, "wordcloud.png", "image/png", key="download_wordcloud")
    else:
        st.write("Brak dostępnych słów kluczowych do wygenerowania chmury słów.")



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
