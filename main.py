import logging
import streamlit as st
from Bio import Entrez
import time
import os
from urllib.error import HTTPError
import pandas as pd
import plotly.express as px
import re
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import io


Entrez.email = "martynapradela@gmail.com"

# Funkcja konfigurująca logger
def setup_logger(term):
    log_folder = f"log-{term}"
    os.makedirs(log_folder, exist_ok=True)
    log_file = os.path.join(log_folder, "search.log")

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ],
    )
    return logging.getLogger(__name__)


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
        
def extract_github_links(abstract):
    """
    Wyszukuje linki GitHub w abstrakcie i usuwa niepożądane znaki na końcu.
    """
    if not abstract:
        return "No GitHub links"

    # Wyrażenie regularne dla linków GitHub (zatrzymuje się na pierwszym znaku kończącym link)
    links = re.findall(r'https://github\.com/[^\s\)]+', abstract)

    # Oczyszczenie linków z potencjalnych nawiasów lub przecinków na końcu
    cleaned_links = [link.rstrip("),.") for link in links]

    return ", ".join(cleaned_links) if cleaned_links else "No GitHub links"       

def fetch_articles_with_details(pmids, start_year, end_year, logger=None):
    """
    Pobiera szczegółowe dane artykułów z abstraktami i informacją o linkach GitHub.
    Wyniki są filtrowane na podstawie zakresu lat publikacji.
    """
    all_articles_data = []
    batch_size = 1000  # Liczba PMIDs przetwarzanych w jednej iteracji efetch
    total_batches = len(pmids) // batch_size + (1 if len(pmids) % batch_size > 0 else 0)  # Liczba partii

    logger.info(f"Rozpoczynanie pobierania szczegółowych danych dla {len(pmids)} PMIDs...")
    for i, start in enumerate(range(0, len(pmids), batch_size)):
        batch_pmids = pmids[start:start + batch_size]
        try:
            # Pobieranie danych w partiach
            stream = Entrez.efetch(db="pubmed", id=",".join(batch_pmids), retmode="xml")
            results = Entrez.read(stream)
            stream.close()

            for article in results.get("PubmedArticle", []):
                try:
                    # Podstawowe informacje
                    pmid = article["MedlineCitation"]["PMID"]
                    title = article["MedlineCitation"]["Article"]["ArticleTitle"]
                    country = article["MedlineCitation"]["MedlineJournalInfo"].get("Country", "Country not found")
                    
                    # Data publikacji
                    pub_date = article["MedlineCitation"]["Article"].get("Journal", {}).get("JournalIssue", {}).get("PubDate", {})
                    pub_date_year = pub_date.get("Year", None)

                    # Filtruj po roku publikacji
                    if pub_date_year and (start_year <= int(pub_date_year) <= end_year):
                        # Autorzy
                        authors = ", ".join([
                            f"{author.get('LastName', '')} {author.get('ForeName', '')}".strip()
                            for author in article["MedlineCitation"]["Article"].get("AuthorList", [])
                        ])

                        # DOI URL
                        e_location_id = article["PubmedData"]["ArticleIdList"]
                        doi_url = "DOI not found"
                        for eloc in e_location_id:
                            if eloc.attributes.get("IdType") == "doi":
                                doi_url = f"https://doi.org/{eloc}"
                                break

                        # Abstrakt i status abstraktu
                        abstract_text = ""
                        has_abstract = "No"
                        if "Abstract" in article["MedlineCitation"]["Article"]:
                            abstract = article["MedlineCitation"]["Article"]["Abstract"]["AbstractText"]
                            if isinstance(abstract, list):
                                abstract_text = " ".join(abstract)
                            else:
                                abstract_text = abstract
                            has_abstract = "Yes"

                        # Linki GitHub i status
                        github_links = extract_github_links(abstract_text)
                        has_github_links = "Yes" if "https://github.com" in github_links else "No"

                        # Kompletowanie danych artykułu
                        article_data = {
                            "PMID": pmid,
                            "Title": title,
                            "Authors": authors,
                            "Country": country,
                            "Publication Date": pub_date_year,
                            "Abstract": abstract_text,
                            "Has Abstract": has_abstract,
                            "Has GitHub Links": has_github_links,
                            "GitHub Links": github_links,
                            "DOI URL": doi_url,
                        }
                        all_articles_data.append(article_data)

                except Exception as e:
                    logger.error(f"Błąd przetwarzania artykułu: {e}")
                    continue

            logger.info(f"Pobrano szczegóły dla {len(batch_pmids)} PMIDs.")

        except Exception as e:
            logger.error(f"Błąd podczas pobierania partii PMIDs: {e}")
            continue

    logger.info(f"Łącznie pobrano szczegóły dla {len(all_articles_data)} artykułów.")
    return all_articles_data


# Funkcje do pobierania wyników
def get_pmids_with_full_pagination(term, start_year=1900, end_year=2023, logger=None):
    total_pmids = []

    logger.info("### Rozpoczynanie pobierania PMIDs ###")
    for year in range(start_year, end_year + 1):
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
            logger.info(f"Pobrano {len(pmids)} PMIDs dla roku {year}.")

            if len(pmids) == 9999:
                monthly_pmids = get_pmids_for_months(term, year, logger)
                yearly_pmids.extend(monthly_pmids)

            total_pmids.extend(yearly_pmids)
        except Exception as e:
            logger.error(f"Błąd podczas pobierania wyników dla roku {year}: {e}")

    unique_pmids = list(set(total_pmids))
    duplicate_count = len(total_pmids) - len(unique_pmids)

    logger.info("### Podsumowanie ###")
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
    for day in range(1, 32):  # Obsługa dni miesiąca
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
    """
    Tworzy interaktywny wykres słupkowy dla liczby publikacji na przestrzeni lat i umożliwia zapis na komputer.
    
    Args:
        data (dict): Słownik, gdzie klucze to lata, a wartości to liczba publikacji.
        output_file (str): Ścieżka do pliku, w którym ma być zapisany wykres (HTML lub PNG).
        logger (Logger): Logger do logowania (opcjonalne).
    
    Returns:
        None: Wykres otwiera się w przeglądarce.
    """
    try:
        # Tworzenie DataFrame z danych
        df = pd.DataFrame({
            "Year": list(data.keys()),
            "Number of Publications": list(data.values())
        })

        # Tworzenie wykresu słupkowego
        fig = px.bar(
            df,
            x="Year",
            y="Number of Publications",
            title="Liczba publikacji na przestrzeni lat",
            labels={"Year": "Rok", "Number of Publications": "Liczba publikacji"},
            hover_data={"Number of Publications": True},
        )

        # Dodanie interaktywności
        fig.update_traces(marker=dict(line=dict(width=1, color="DarkSlateGrey")))
        fig.update_layout(
            hovermode="x",
            xaxis=dict(tickangle=45),
            template="plotly_white"
        )

        # Wyświetlenie wykresu
        fig.show()
        if logger:
            logger.info("Wykres słupkowy został wygenerowany pomyślnie.")

        # Zapis do pliku
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
    import geopandas as gpd
    from matplotlib.colors import LinearSegmentedColormap

    # Zliczanie publikacji według krajów
    country_counts = dataframe['Country'].value_counts().reset_index()
    country_counts.columns = ['Country', 'Count']

    # Załadowanie kształtów krajów z lokalnego pliku .shp
    world = gpd.read_file("C:/Users/Martyna/OneDrive/Pulpit/web scraping/ne_110m_admin_0_countries.shp")

    # Dopasowanie liczby publikacji do kształtów krajów
    world = world.merge(country_counts, left_on='NAME', right_on='Country', how='left')

    # Wypełnianie braków wartości zerami
    world['Count'] = world['Count'].fillna(0)

    # Tworzenie niestandardowej palety kolorów
    custom_cmap = LinearSegmentedColormap.from_list(
        "custom_cmap",
        ["grey", "yellow", "blue", "green"]
    )

    # Rysowanie mapy
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
    return fig  # Zwracamy obiekt `Figure`


# Streamlit interfejs
st.title("PubMed Data Viewer HEJ")
term = st.text_input("Wprowadź termin wyszukiwania:", value="heart")
start_year = st.number_input("Podaj początkowy rok:", value=2000, step=1)
end_year = st.number_input("Podaj końcowy rok:", value=2023, step=1)

# Inicjalizacja session state
if "articles_data" not in st.session_state:
    st.session_state["articles_data"] = None
if "publications_by_year" not in st.session_state:
    st.session_state["publications_by_year"] = None

if st.button("Wyszukaj"):
    is_data_loaded = False
    with st.spinner("Pobieranie danych..."):
        try:
            logger = setup_logger(term)

            # Pobierz PMIDs z podanego zakresu lat
            unique_pmids, duplicate_count = get_pmids_with_full_pagination(term, start_year, end_year, logger)

            # Pobierz szczegóły artykułów na podstawie PMIDs
            articles_data = fetch_articles_with_details(unique_pmids, start_year, end_year, logger)
            st.session_state["articles_data"] = articles_data

            # Generowanie danych do wykresu słupkowego
            publications_by_year = {}
            for article in articles_data:
                pub_date = article["Publication Date"]
                if pub_date.isdigit():
                    year = int(pub_date)  # Konwersja na liczbę
                    publications_by_year[year] = publications_by_year.get(year, 0) + 1
            st.session_state["publications_by_year"] = publications_by_year

            is_data_loaded = True
        except Exception as e:
            st.error(f"Wystąpił błąd podczas pobierania danych: {e}")
            is_data_loaded = False

# Wyświetlanie danych tylko wtedy, gdy istnieją w session state
if st.session_state["articles_data"]:
    # Wyświetlenie tabeli
    st.write("### Szczegółowe dane artykułów ###")
    df = pd.DataFrame(st.session_state["articles_data"])
    st.dataframe(df)

    # Wyświetlenie wykresu słupkowego
    st.subheader(f"Wykres słupkowy dla: {term} w PubMed")
    publications_by_year = st.session_state["publications_by_year"]

    if publications_by_year:
        fig_bar = px.bar(
            x=list(publications_by_year.keys()),
            y=list(publications_by_year.values()),
            labels={"x": "Rok", "y": "Liczba publikacji"},
            title=f"Wykres słupkowy dla: {term} w PubMed"
        )
        fig_bar.update_traces(hovertemplate='Rok: %{x}<br>Liczba publikacji: %{y}<extra></extra>')
        st.plotly_chart(fig_bar, use_container_width=True)

        # Zapisanie wykresu słupkowego jako PNG
        buf_bar = io.BytesIO()
        fig_bar.write_image(buf_bar, format="png", scale=2)  # Używamy Plotly do zapisania jako PNG
        buf_bar.seek(0)

        # Pobranie wykresu słupkowego jako PNG
        st.download_button(
            label="Pobierz wykres słupkowy jako PNG",
            data=buf_bar,
            file_name="wykres_slupkowy.png",
            mime="image/png",
            key="download_bar_chart"
        )

    # Wyświetlenie mapy świata
    st.write("### Mapa publikacji według krajów ###")
    map_figure = create_map(df)
    st.pyplot(map_figure)

    # Zapisanie mapy jako PNG
    buf_map = io.BytesIO()
    map_figure.savefig(buf_map, format="png", dpi=300, bbox_inches="tight")
    buf_map.seek(0)

    # Pobranie mapy świata jako PNG
    st.download_button(
        label="Pobierz mapę jako PNG",
        data=buf_map,
        file_name="mapa_publikacji.png",
        mime="image/png",
        key="download_map_chart"
    )






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
