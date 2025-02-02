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
from bs4 import BeautifulSoup
from urllib.parse import unquote


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
        keywords_list = article["MedlineCitation"].get("KeywordList", [])

        clean_keywords = []
        for sublist in keywords_list:
            for kw in sublist:
                clean_text = BeautifulSoup(str(kw), "html.parser").get_text()  # Usunięcie tagów HTML
                clean_keywords.append(clean_text.strip())

        return "; ".join(clean_keywords) if clean_keywords else "Brak słów kluczowych"
    except Exception as e:
        return "Błąd pobierania"

def generate_wordcloud(text):
    if text and text != "Brak słów kluczowych":
        wordcloud = WordCloud(width=800, height=400, background_color="white").generate(text)

        fig, ax = plt.subplots(figsize=(10, 5))
        ax.imshow(wordcloud, interpolation="bilinear")
        ax.axis("off")
        return fig
    return None 
        
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

                        # Pobieranie słów kluczowych
                        keywords = set()
                        if "KeywordList" in article["MedlineCitation"]:
                            for kw_list in article["MedlineCitation"]["KeywordList"]:
                                for kw in kw_list:
                                    clean_text = BeautifulSoup(str(kw), "html.parser").get_text()
                                    clean_text = unquote(clean_text)  
                                    keywords.add(clean_text.strip())  

                        # Pobieranie MeshHeadingList
                        if "MeshHeadingList" in article["MedlineCitation"]:
                            mesh_terms = article["MedlineCitation"]["MeshHeadingList"]
                            keywords.update(str(term["DescriptorName"]) for term in mesh_terms)

                        all_keywords.update(keywords)  
                        keywords_text = "; ".join(keywords) if keywords else "No Keywords"

                        # Pobieranie abstraktu i linków GitHub
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

# Nagłówek aplikacji
st.title("PubMed Data Viewer")

# Formularz wyszukiwania
term = st.text_input("Wprowadź termin wyszukiwania:", value="heart cancer", key="search_term_input")
start_year = st.number_input("Podaj początkowy rok:", value=2000, step=1, key="start_year_input")
end_year = st.number_input("Podaj końcowy rok:", value=2023, step=1, key="end_year_input")

# Inicjalizacja zmiennych w session_state
if "articles_data" not in st.session_state:
    st.session_state["articles_data"] = None
if "publications_by_year" not in st.session_state:
    st.session_state["publications_by_year"] = None
if "github_repos_by_year" not in st.session_state:
    st.session_state["github_repos_by_year"] = None
if "all_keywords" not in st.session_state:
    st.session_state["all_keywords"] = None
if "wordcloud_fig" not in st.session_state:
    st.session_state["wordcloud_fig"] = None  

# Pobieranie danych po kliknięciu przycisku
if st.button("Wyszukaj"):
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
            all_keywords = []

            for article in articles_data:
                pub_date = article["Publication Date"]
                has_github = article["Has GitHub Links"] == "Yes"
                keywords = article["Keywords"].split("; ") if article["Keywords"] != "No Keywords" else []
                all_keywords.extend(keywords)

                if pub_date.isdigit():
                    publications_by_year[pub_date] = publications_by_year.get(pub_date, 0) + 1
                    if has_github:
                        github_repos_by_year[pub_date] = github_repos_by_year.get(pub_date, 0) + 1

            st.session_state["articles_data"] = articles_data
            st.session_state["publications_by_year"] = publications_by_year
            st.session_state["github_repos_by_year"] = github_repos_by_year

            if not st.session_state["all_keywords"]:
                st.session_state["all_keywords"] = " ".join(all_keywords) if all_keywords else "No Keywords"

            if st.session_state["all_keywords"] and st.session_state["all_keywords"] != "No Keywords":
                wordcloud = WordCloud(width=800, height=400, background_color="white").generate(st.session_state["all_keywords"])
                fig, ax = plt.subplots(figsize=(10, 5))
                ax.imshow(wordcloud, interpolation="bilinear")
                ax.axis("off")
                st.session_state["wordcloud_fig"] = fig
            else:
                st.session_state["wordcloud_fig"] = None  

        except Exception as e:
            st.error(f"Wystąpił błąd: {e}")

# Jeśli są dane, wyświetl wykresy
if st.session_state["articles_data"]:
    df = pd.DataFrame(st.session_state["articles_data"])
    unique_countries = sorted(df["Country"].dropna().unique())
    unique_countries.insert(0, "Wszystkie kraje")
    all_keywords_list = sorted(set("; ".join(df["Keywords"].dropna()).split("; ")))

    ## Filtracja po kraju**
    st.subheader("Filtracja po kraju")
    selected_country = st.selectbox("Wybierz kraj:", unique_countries, index=0, key="country_filter")

    df_filtered = df if selected_country == "Wszystkie kraje" else df[df["Country"] == selected_country]
    chart_title = "Liczba publikacji na przestrzeni lat" if selected_country == "Wszystkie kraje" else f"Liczba publikacji dla {selected_country}"
    publications_by_year = df_filtered["Publication Date"].value_counts().sort_index()

    if not publications_by_year.empty:
        fig_country = px.bar(x=publications_by_year.index.astype(int), y=publications_by_year.values,
                             labels={"x": "Rok", "y": "Liczba publikacji"}, title=chart_title,
                             color_discrete_sequence=["deepskyblue"])
        fig_country.update_xaxes(type="linear", dtick=1)
        st.plotly_chart(fig_country, use_container_width=True)

        buf_country = io.BytesIO()
        fig_country.write_image(buf_country, format="png", scale=2)
        buf_country.seek(0)
        st.download_button("Pobierz wykres jako PNG", buf_country, f"wykres_kraj.png", "image/png", key="download_country")

    #Filtracja po słowach kluczowych**
    st.subheader("Filtracja po słowach kluczowych")
    selected_keywords = st.multiselect("Wybierz słowa kluczowe:", all_keywords_list, key="keyword_filter")

    df_keyword_filtered = df[df["Keywords"].apply(lambda x: any(kw in x for kw in selected_keywords))] if selected_keywords else df
    publications_by_year_keyword = df_keyword_filtered["Publication Date"].value_counts().sort_index()

    keyword_title = f"Liczba publikacji według słów kluczowych: {', '.join(selected_keywords)}" if selected_keywords else "Liczba publikacji według słów kluczowych"

    if not publications_by_year_keyword.empty:
        fig_keyword = px.bar(
            x=publications_by_year_keyword.index.astype(int), 
            y=publications_by_year_keyword.values,
            labels={"x": "Rok", "y": "Liczba publikacji"}, 
            title=keyword_title,
            color_discrete_sequence=["deepskyblue"]
        )
        st.plotly_chart(fig_keyword, use_container_width=True)

        
        buf_keyword = io.BytesIO()
        fig_keyword.write_image(buf_keyword, format="png", scale=2)
        buf_keyword.seek(0)
        st.download_button("Pobierz wykres jako PNG", buf_keyword, "filtracja_slow_kluczowych.png", "image/png", key="download_keyword_chart")


    ##Filtracja po kraju i słowach kluczowyc
    st.subheader("Filtracja po kraju i słowach kluczowych")
    selected_country_2 = st.selectbox("Wybierz kraj:", unique_countries, index=0, key="country_keyword_filter")
    selected_keywords_2 = st.multiselect("Wybierz słowa kluczowe:", all_keywords_list, key="keyword_country_filter")

    df_combined_filtered = df if selected_country_2 == "Wszystkie kraje" else df[df["Country"] == selected_country_2]
    if selected_keywords_2:
        df_combined_filtered = df_combined_filtered[df_combined_filtered["Keywords"].apply(lambda x: any(kw in x for kw in selected_keywords_2))]

    publications_by_year_combined = df_combined_filtered["Publication Date"].value_counts().sort_index()

    combined_title = f"Liczba publikacji według kraju: {selected_country_2} i słów kluczowych: {', '.join(selected_keywords_2)}" if selected_keywords_2 else f"Liczba publikacji według kraju: {selected_country_2}"

    if not publications_by_year_combined.empty:
        fig_combined = px.bar(
            x=publications_by_year_combined.index.astype(int), 
            y=publications_by_year_combined.values,
            labels={"x": "Rok", "y": "Liczba publikacji"}, 
            title=combined_title,
            color_discrete_sequence=["deepskyblue"]
        )
        st.plotly_chart(fig_combined, use_container_width=True)

        
        buf_combined = io.BytesIO()
        fig_combined.write_image(buf_combined, format="png", scale=2)
        buf_combined.seek(0)
        st.download_button("Pobierz wykres jako PNG", buf_combined, "filtracja_kraj_slow_kluczowych.png", "image/png", key="download_combined_chart")

    #Wykres repozytoriów GitHub
    st.subheader("Liczba znalezionych repozytoriów GitHub w określonym przedziale czasowym")
    github_repos_by_year = st.session_state["github_repos_by_year"]

    if github_repos_by_year:
        fig_github = px.bar(x=list(map(int, github_repos_by_year.keys())), y=list(github_repos_by_year.values()),
                            labels={"x": "Rok", "y": "Liczba repozytoriów GitHub"}, title="Liczba repozytoriów GitHub w określonym przedziale czasowym",
                            color_discrete_sequence=["deepskyblue"])
        fig_github.update_xaxes(type="linear", dtick=1)
        st.plotly_chart(fig_github, use_container_width=True)

        buf_github = io.BytesIO()
        fig_github.write_image(buf_github, format="png", scale=2)
        buf_github.seek(0)
        st.download_button("Pobierz wykres jako PNG", buf_github, f"wykres_github.png", "image/png", key="download_github")


    # Mapa świata 
    st.write("### Mapa globalnego udziału krajów w publikacjach naukowych ###")
    map_figure = create_map(df)
    st.pyplot(map_figure)

    buf_map = io.BytesIO()
    map_figure.savefig(buf_map, format="png", dpi=300, bbox_inches="tight")
    buf_map.seek(0)
    st.download_button("Pobierz mapę jako PNG", buf_map, "mapa_publikacji.png", "image/png", key="download_map_chart")


    #TOP 5
    st.write(f"### TOP 5 krajów z największym udziałem w tworzeniu zapytania: {term}")
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

    st.write("### Chmura słów kluczowych ###")

    if st.session_state["wordcloud_fig"]:
        st.pyplot(st.session_state["wordcloud_fig"])  

        buf_wordcloud = io.BytesIO()
        st.session_state["wordcloud_fig"].savefig(buf_wordcloud, format="png", dpi=300, bbox_inches="tight")
        buf_wordcloud.seek(0)
        st.download_button("Pobierz chmurę słów jako PNG", buf_wordcloud, "wordcloud.png", "image/png")
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
