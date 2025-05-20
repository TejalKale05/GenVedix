# GenVedix - An Aisle to Gene World
import streamlit as st
import pandas as pd
import requests
import base64
import io 
import os
import time
from Bio import Entrez, SeqIO 
from PIL import Image 
import networkx as nx
from pyvis.network import Network
import json
from Bio.Seq import UndefinedSequenceError, Seq
from Bio.Data import CodonTable
import plotly.graph_objects as go
import plotly.express as px
from collections import Counter
from scipy.stats import zscore

# --- Primer3 import for Primer Design ---
try:
    import primer3
except ImportError:
    primer3 = None
    st.sidebar.warning("`primer3-py` library not found. Primer Design feature will not work. Install: `pip install primer3-py`")

# --- py3Dmol import for 3D Structure Visualization ---
try:
    import py3Dmol
except ImportError:
    py3Dmol = None
    st.sidebar.warning("`py3Dmol` library not found. 3D Structure Viewer will not work. Install: `pip install py3Dmol`")


# --- Imports for New Features ---
from Bio import Phylo


# --- SIMPLE, ELEGANT, HIGH-CONTRAST CSS ---
def get_base64_of_bin_file(bin_file):
    with open(bin_file, 'rb') as f:
        data = f.read()
    return base64.b64encode(data).decode()


st.markdown(f'''<style>
body, .stApp {{
    background: #21618C;
    color: #eaffea;
    font-family: 'Segoe UI', 'Roboto', 'Arial', sans-serif;
    font-size: 1.18em;
}}
.stApp {{
    background-image: none !important;
}}
.stSidebar {{
    background: #4CAF50 !important;
    color: #0a3a4a !important;
    border-right: 2px solid #4CAF50;
    font-size: 1.15em !important;
}}
[data-testid="stSidebar"] .css-1v0mbdj,
[data-testid="stSidebar"] .css-1d391kg,
[data-testid="stSidebar"] label,
[data-testid="stSidebar"] .stRadio label,
[data-testid="stSidebar"] .stSelectbox label,
[data-testid="stSidebar"] .stTextInput label,
[data-testid="stSidebar"] .stExpanderHeader,
[data-testid="stSidebar"] .st-bb,
[data-testid="stSidebar"] .st-c3,
[data-testid="stSidebar"] .st-c4 {{
    color: #0a3a4a !important;
    font-weight: 600;
    font-size: 1.13em !important;
}}
[data-testid="stSidebar"] .css-10trblm,
[data-testid="stSidebar"] h1,
[data-testid="stSidebar"] h2 {{
    color: #0a3a4a !important;
    text-shadow: none !important;
    font-size: 2.1em !important;
}}
[data-testid="stTabs"] button {{
    color: #fff !important;
}}
[data-testid="stTabs"] button[aria-selected="true"] {{
    color: #fff !important;
    border-bottom: 3px solid #fff !important;
}}
[data-testid="stTabs"] svg {{
    color: #fff !important;
    fill: #fff !important;
}}
.stButton>button, .stDownloadButton>button, .custom-link-button {{
    background: linear-gradient(90deg, #4CAF50 0%, #0a3a4a 100%) !important;
    color: #fff !important;
    border-radius: 8px;
    border: none;
    font-weight: 600;
    padding: 0.6em 1.7em;
    margin: 0.5em 0;
    box-shadow: 0 2px 8px rgba(76,175,80,0.15);
    transition: background 0.3s, box-shadow 0.3s;
    font-size: 1.13em !important;
}}
.stButton>button:hover, .stDownloadButton>button:hover, .custom-link-button:hover {{
    background: linear-gradient(90deg, #0a3a4a 0%, #4CAF50 100%) !important;
    box-shadow: 0 4px 16px rgba(0,188,212,0.18);
}}
h1, .main-header {{
    font-family: 'Segoe UI', 'Roboto', Arial, sans-serif;
    font-size: 2.7rem;
    letter-spacing: 1.2px;
    color: #4CAF50;
    text-shadow: 0 2px 8px #0a3a4a;
}}
h2, h3, h4 {{ /* Applied to general h2,h3,h4 if not overridden by specific classes */
    color: #0a3a4a;
    font-family: 'Segoe UI', 'Roboto', Arial, sans-serif;
    letter-spacing: 1px;
    font-size: 1.35em !important;
}}
/* Specific styling for Home Page section headers */
.home-section-header {{
    color: #4CAF50 !important;
    font-size: 2.2em !important;
    border-bottom: 2.5px solid #00bcd4;
    padding-bottom: 0.4em;
    margin-top: 1.5em;
    margin-bottom: 1em;
    text-align: center;
}}
/* Styling for expander headers on the home page for tool details */
.home-tool-expander .stExpanderHeader {{
    color: #0a3a4a !important; /* Main accent color for expander title */
    font-weight: 700 !important;
    font-size: 1.25em !important; /* Slightly larger for home page context */
    background-color: #2a2e35; /* Darker background for the header itself */
    border-radius: 8px 8px 0 0; /* Rounded top corners */
    padding: 0.8em 1em !important;
    border-bottom: 1.5px solid #4CAF50; /* Accent border */
}}
.home-tool-expander > div > div:nth-child(2) {{ /* Targeting the content part of the expander */
    background-color: #23272f; /* Main card background for content */
    border: 1.5px solid #4CAF50;
    border-top: none; /* Remove top border as header has bottom border */
    border-radius: 0 0 8px 8px; /* Rounded bottom corners */
    padding: 1.5em !important;
    font-size: 1.05em;
    color: #e0e0e0;
}}
.home-tool-expander .stExpanderHeader p {{ /* Ensure text inside expander header is styled */
    color: #00bcd4 !important;
}}
hr {{
    border-top: 1.5px solid #00bcd4;
}}
.stExpanderHeader {{
    color: #00bcd4 !important;
    font-weight: 700;
    font-size: 1.13em !important;
}}
.result-box, .feature-card, .main-container {{
    background: #23272f;
    border-radius: 12px;
    border: 1.5px solid #4CAF50;
    box-shadow: 0 2px 12px rgba(0,188,212,0.08);
    padding: 1.5em 1.5em 1em 1.5em;
    margin-bottom: 1.5em;
    font-size: 1.13em;
}}
.stDownloadButton>button, .stButton>button {{
    font-size: 1.13em;
}}
.stSelectbox>div>div {{
    background: #4CAF50 !important;
    color: #0a3a4a !important;
    font-size: 1.13em !important;
}}
.stTextArea>div>textarea {{
    background: #4CAF50 !important;
    color: #0a3a4a !important;
    border: 1.5px solid #00bcd4 !important;
    border-radius: 6px;
    font-size: 1.13em !important;
}}
.stAlert, .stWarning, .stError, .stSuccess {{
    border-radius: 8px !important;
    font-size: 1.13em !important;
}}
::-webkit-scrollbar {{
    width: 10px;
    background: #23272f;
}}
::-webkit-scrollbar-thumb {{
    background: #4CAF50;
    border-radius: 8px;
}}
label, .stTextInput label, .stSelectbox label, .stRadio label {{
    color: #4CAF50 !important;
    font-weight: 600;
    font-size: 1.18em !important;
}}
input[type="text"], input[type="search"], input[type="email"], input[type="password"], textarea, .stTextInput>div>input {{
    background: #4CAF50 !important;
    color: #0a3a4a !important;
    border: 1.5px solid #00bcd4 !important;
    border-radius: 6px !important;
    font-size: 1.13em !important;
}}
.tool-card-grid-container {{
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(300px, 1fr));
    gap: 25px;
    padding: 15px 0;
}}
.tool-card-style {{
    background-color: #2a2e35;
    border: 1.5px solid #00bcd4;
    border-radius: 12px;
    padding: 1.5em;
    display: flex;
    flex-direction: column;
    align-items: center;
    text-align: center;
    box-sizing: border-box;
    min-height: 300px;
    justify-content: space-between;
    box-shadow: 0 2px 10px rgba(0,188,212,0.1);
    height: 100%;
}}
.tool-card-title {{
    color: #4CAF50 !important;
    margin-top: 0.2em;
    margin-bottom: 0.6em;
    font-size: 1.4em !important;
    line-height: 1.3;
    font-weight: 700;
}}
.tool-card-description {{
    font-size: 0.9em;
    color: #b0bec5;
    margin-bottom: 1em;
    flex-grow: 1;
    line-height: 1.5;
    text-align: left;
    width: 100%;
}}
.tool-card-button-wrapper {{
    width: 100%;
    display: flex;
    justify-content: center;
    margin-top: auto;
}}
.tool-card-button-wrapper .stButton>button {{
    background: linear-gradient(90deg, #4CAF50 0%, #00bcd4 100%) !important;
    color: #fff !important;
    border-radius: 8px !important;
    border: none !important;
    font-weight: 600 !important;
    padding: 0.6em 1.5em !important;
    box-shadow: 0 2px 8px rgba(76,175,80,0.15) !important;
    transition: background 0.3s, box-shadow 0.3s !important;
    font-size: 1.0em !important;
    width: clamp(160px, 85%, 230px) !important;
}}
.tool-card-button-wrapper .stButton>button:hover {{
    background: linear-gradient(90deg, #00bcd4 0%, #4CAF50 100%) !important;
    box-shadow: 0 4px 16px rgba(0,188,212,0.18) !important;
}}
.section-subheader {{ /* General section subheader for other pages */
    color: #4CAF50;
    font-size: 1.9em !important;
    border-bottom: 2.5px solid #00bcd4;
    padding-bottom: 0.4em;
    margin-top: 2em;
    margin-bottom: 1.2em;
}}

/* Styling for lists within home page expanders */
.home-tool-expander ul {{
    list-style-type: disc;
    padding-left: 25px;
    margin-top: 0.5em;
    margin-bottom: 1em;
    color: #c0c5ce; /* Lighter text for list items */
}}
.home-tool-expander li {{
    margin-bottom: 0.6em;
    line-height: 1.6;
}}
.home-tool-expander strong {{
    color: #4CAF50; /* Accent color for strong/bold text */
}}
.home-tool-expander code {{ /* Styling for inline code/PDB IDs etc. */
    background-color: #333942;
    color: #00bcd4;
    padding: 0.2em 0.4em;
    border-radius: 4px;
    font-family: 'Courier New', Courier, monospace;
    font-size: 0.95em;
}}

</style>''', unsafe_allow_html=True)

st.sidebar.markdown("""
    <h1 style='color:#0a3a4a; font-size:3.2em; font-weight:900; letter-spacing:2px; margin-bottom:0.12em;'>GenVedix</h1>
    <div style='color:#0a3a4a; font-size:1.1em; font-weight:800; margin-bottom:0.7em; margin-top:-0.35em; white-space:nowrap;'>An Aisle To Gene World</div>
""", unsafe_allow_html=True)

ORGANISM_TO_TAXID = {
    "Homo sapiens": 9606, "Mus musculus": 10090, "Drosophila melanogaster": 7227,
    "Caenorhabditis elegans": 6239, "Saccharomyces cerevisiae": 4932,
    "Danio rerio": 7955, "Arabidopsis thaliana": 3702,
    "Rattus norvegicus": 10116, "Gallus gallus": 9031
}
ORGANISM_TO_KEGG_CODE = {
    "Homo sapiens": "hsa", "Mus musculus": "mmu", "Drosophila melanogaster": "dme",
    "Caenorhabditis elegans": "cel", "Saccharomyces cerevisiae": "sce",
    "Danio rerio": "dre", "Arabidopsis thaliana": "ath",
    "Rattus norvegicus": "rno", "Gallus gallus": "gga"
}
ORGANISM_TO_ENSEMBL_SPECIES = {
    "Homo sapiens": "homo_sapiens", "Mus musculus": "mus_musculus",
    "Drosophila melanogaster": "drosophila_melanogaster", "Caenorhabditis elegans": "caenorhabditis_elegans",
    "Saccharomyces cerevisiae": "saccharomyces_cerevisiae", "Danio rerio": "danio_rerio",
    "Arabidopsis thaliana": "arabidopsis_thaliana",
    "Rattus norvegicus": "rattus_norvegicus", "Gallus gallus": "gallus_gallus"
}
ENSEMBL_SERVER = "https://rest.ensembl.org"
ORGANISM_OPTIONS = list(ORGANISM_TO_TAXID.keys())

@st.cache_data(ttl=3600)
def fetch_kegg_pathway_image(pathway_id):
    img_url = f"https://rest.kegg.jp/get/{pathway_id}/image"
    response = requests.get(img_url)
    if response.status_code == 200:
        return response.content, None
    else:
        return None, f"Could not fetch pathway image for {pathway_id} (Status: {response.status_code})"

@st.cache_data(ttl=3600)
def fetch_kegg_pathways(gene_symbol, organism_name):
    try:
        if organism_name not in ORGANISM_TO_KEGG_CODE:
            return None, f"Organism '{organism_name}' not supported for KEGG."
        kegg_code = ORGANISM_TO_KEGG_CODE[organism_name]
        gene_query = f"{kegg_code}:{gene_symbol}"
        pathway_url_direct = f"https://rest.kegg.jp/link/pathway/{gene_query}"
        response_direct = requests.get(pathway_url_direct)
        gene_id_to_use = None
        if response_direct.status_code == 200 and response_direct.text.strip():
            gene_id_to_use = gene_query
        else:
            gene_list_url = f"https://rest.kegg.jp/list/{kegg_code}"
            response_list = requests.get(gene_list_url)
            if response_list.status_code != 200:
                return None, f"Error fetching KEGG gene list for {kegg_code}: {response_list.text}"
            found_gene_id = None
            for line in response_list.text.strip().split('\n'):
                if line:
                    parts = line.split('\t', 1)
                    if len(parts) == 2:
                        current_kegg_gene_id, description = parts
                        desc_parts = description.split(';')
                        gene_symbols_in_desc = desc_parts[0].strip()
                        normalized_gene_symbol = gene_symbol.lower()
                        symbols_to_check = [s.strip().lower() for s in gene_symbols_in_desc.split(',')]
                        if normalized_gene_symbol in symbols_to_check or \
                           f"({normalized_gene_symbol})" in gene_symbols_in_desc.lower() or \
                           f"; {normalized_gene_symbol};" in f"; {gene_symbols_in_desc.lower()};":
                            found_gene_id = current_kegg_gene_id
                            break
            if not found_gene_id:
                 for line in response_list.text.strip().split('\n'):
                    if line:
                        parts = line.split('\t', 1)
                        if len(parts) == 2:
                            current_kegg_gene_id, description = parts
                            if gene_symbol.lower() in description.lower():
                                found_gene_id = current_kegg_gene_id
                                break
            if not found_gene_id:
                return None, f"Gene '{gene_symbol}' not found in KEGG database for organism '{organism_name}'. (Tried exact and broader description search)"
            gene_id_to_use = found_gene_id
        pathway_url = f"https://rest.kegg.jp/link/pathway/{gene_id_to_use}"
        response_pathways = requests.get(pathway_url)
        if response_pathways.status_code != 200:
            return None, f"Error fetching KEGG pathways for {gene_id_to_use}: {response_pathways.text}"
        pathways = []
        for line in response_pathways.text.strip().split('\n'):
            if line:
                parts = line.split('\t', 1)
                if len(parts) == 2:
                    _, pathway_id_full = parts
                    pathway_id = pathway_id_full.replace('path:', '')
                    details_url = f"https://rest.kegg.jp/get/{pathway_id}"
                    details_response = requests.get(details_url)
                    name = "N/A"
                    if details_response.status_code == 200:
                        details_text = details_response.text
                        for detail_line in details_text.split('\n'):
                            if detail_line.startswith('NAME'):
                                name = detail_line.split('NAME', 1)[1].strip().split(' - ')[0]
                                break
                    pathways.append({"Pathway ID": pathway_id, "Name": name, "URL": f"https://www.kegg.jp/pathway/{pathway_id}"})
        if not pathways:
            return None, f"No KEGG pathways found for gene '{gene_symbol}' (KEGG ID: {gene_id_to_use})."
        return pathways, None
    except Exception as e:
        return None, f"Unexpected error in KEGG pathway fetch: {str(e)}"

# --- Reactome Pathways ---
@st.cache_data(ttl=3600)
def fetch_reactome_pathways(gene_symbol, organism_name):
    try:
        # Try Ensembl ID as well as symbol
        species = ORGANISM_TO_ENSEMBL_SPECIES.get(organism_name, None)
        ids_to_try = [gene_symbol]
        if species:
            url = f"https://rest.ensembl.org/lookup/symbol/{species}/{gene_symbol}"
            headers = {"Content-Type": "application/json"}
            r = requests.get(url, headers=headers)
            if r.status_code == 200:
                ensembl_id = r.json().get("id")
                if ensembl_id:
                    ids_to_try.append(ensembl_id)
        for query_id in ids_to_try:
            mapping_url = f"https://reactome.org/ContentService/data/mapping/ENSEMBL/{query_id}/pathways"
            r = requests.get(mapping_url)
            if r.status_code == 200 and isinstance(r.json(), list) and len(r.json()) > 0:
                pathways = []
                for entry in r.json():
                    pathways.append({
                        "Pathway ID": entry.get("stId", ""),
                        "Name": entry.get("displayName", ""),
                        "URL": f"https://reactome.org/content/detail/{entry.get('stId', '')}"
                    })
                return pathways, None
        return None, f"No Reactome pathways found for gene '{gene_symbol}' (tried symbol and Ensembl ID)."
    except Exception as e:
        return None, f"Error searching Reactome: {str(e)}"

@st.cache_data(ttl=3600*6)
def fetch_gene_annotation(gene_name, organism):
    Entrez.email = "kaletejal05@mail.com"
    search_term = f"{gene_name}[Gene Name] AND {organism}[Organism]"
    try:
        handle = Entrez.esearch(db="gene", term=search_term, retmax=5)
        record = Entrez.read(handle)
        handle.close()
    except Exception as e:
        return [], f"NCBI ESearch failed: {e}"
    results = []
    if not record["IdList"]:
        return results, f"No gene ID found for '{gene_name}' in '{organism}'."
    for gene_id in record["IdList"]:
        try:
            summary_handle = Entrez.esummary(db="gene", id=gene_id)
            summary_record = Entrez.read(summary_handle)
            summary_handle.close()
        except Exception as e:
            st.warning(f"NCBI ESummary failed for Gene ID {gene_id}: {e}")
            continue
        if not summary_record['DocumentSummarySet']['DocumentSummary']: continue
        docsum = summary_record['DocumentSummarySet']['DocumentSummary'][0]
        nucleotide_ids = "N/A"
        if 'GenomicInfo' in docsum and docsum['GenomicInfo']:
            genomic_info_list = docsum['GenomicInfo']
            if isinstance(genomic_info_list, list) and genomic_info_list:
                acc_versions = [info.get('ChrAccVer', None) for info in genomic_info_list if 'ChrAccVer' in info]
                acc_versions = [acc for acc in acc_versions if acc]
                if acc_versions: nucleotide_ids = ", ".join(list(set(acc_versions)))
                else:
                    other_ids = [info.get('AssemblyAcc', None) for info in genomic_info_list if 'AssemblyAcc' in info]
                    other_ids = [oid for oid in other_ids if oid]
                    if other_ids: nucleotide_ids = ", ".join(list(set(other_ids)))
        elif 'OtherNames' in docsum:
            other_names_str = docsum['OtherNames']
            found_nuc_ids = [name for name in other_names_str.split(',') if name.strip().startswith(("NM_", "NC_", "XM_", "XR_"))]
            if found_nuc_ids: nucleotide_ids = ", ".join(list(set(id.strip() for id in found_nuc_ids)))
        results.append({
            "Gene Symbol": docsum.get('Name', gene_name), "Official Name": docsum.get('Description', 'N/A'),
            "Function": docsum.get('Summary', 'No function description available'), "Chromosome": docsum.get('Chromosome', 'N/A'),
            "Position": docsum.get('MapLocation', 'N/A'), "Exons": docsum.get('ExonCount', 'N/A'),
            "Source Organism": organism, "Gene ID": gene_id,
            "Annotation Status": "Complete" if docsum.get('Status', '').lower() in ["reviewed", "validated"] else "Partial",
            "Aliases": docsum.get('OtherAliases', 'N/A'), "Nucleotide IDs": nucleotide_ids
        })
    if not results: return [], f"Could not retrieve detailed summaries for gene IDs found for '{gene_name}'."
    return results, None

@st.cache_data(ttl=3600*6)
def fetch_gene_sequence(nuccore_id):
    Entrez.email = "kaletejal05@mail.com"
    if not nuccore_id or nuccore_id == "N/A": return None, "No valid Nucleotide ID provided."
    try:
        handle = Entrez.efetch(db="nucleotide", id=nuccore_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
    except Exception as e: return None, f"Error fetching GenBank record for {nuccore_id}: {e}"
    features = []
    for feat_idx, feat in enumerate(record.features):
        if feat_idx >= 10: break
        feature_info = {"type": feat.type, "location": str(feat.location)}
        if feat.qualifiers:
            quals = {k: v[0] if len(v)==1 else v for k,v in feat.qualifiers.items()}
            feature_info["qualifiers"] = str(quals)[:200] + "..." if len(str(quals)) > 200 else str(quals)
        features.append(feature_info)
    try:
        seq_str = str(record.seq)[:1000]; seq_len = len(record.seq)
    except UndefinedSequenceError: seq_str = "[Sequence content is large or constructed...]"; seq_len = 0
    except AttributeError: seq_str = "[Sequence not available...]"; seq_len = 0
    return {"id": record.id, "description": record.description, "length": seq_len, "sequence": seq_str, "features": features}, None

@st.cache_data(ttl=3600*6)
def fetch_fasta_sequence_data(nuccore_id):
    Entrez.email = "kaletejal05@mail.com"
    if not nuccore_id or nuccore_id == "N/A": return None, "No valid Nucleotide ID provided."
    try:
        handle = Entrez.efetch(db="nucleotide", id=nuccore_id, rettype="fasta", retmode="text")
        fasta_str = handle.read()
        handle.close()
        return fasta_str, None
    except Exception as e: return None, f"Error fetching FASTA for {nuccore_id}: {e}"

@st.cache_data(ttl=3600*2)
def fetch_string_interaction_data_for_pyvis(gene_symbol, organism, score_threshold, limit):
    species_id = ORGANISM_TO_TAXID.get(organism)
    if not species_id: return pd.DataFrame(), f"Organism '{organism}' not supported by STRING DB via this app."
    api_url = "https://string-db.org/api"; output_format = "tsv-no-header"; method = "network"
    request_url = f"{api_url}/{output_format}/{method}"
    params = {"identifiers": gene_symbol.upper(), "species": species_id, "required_score": score_threshold, "limit": limit, "caller_identity": "GenVedix_App_v1"}
    try:
        response = requests.get(request_url, params=params, timeout=20)
        response.raise_for_status()
        if not response.text.strip(): return pd.DataFrame(columns=['protein1', 'protein2', 'score']), f"No interactions found for {gene_symbol} with score >= {score_threshold}."
        data = io.StringIO(response.text)
        interactions_df = pd.read_csv(data, sep='\t', usecols=[2, 3, 5], names=['protein1', 'protein2', 'score'], header=None)
        if interactions_df.empty: return pd.DataFrame(columns=['protein1', 'protein2', 'score']), f"No interactions parsed for {gene_symbol}."
        return interactions_df, None
    except requests.exceptions.HTTPError as e: return pd.DataFrame(), f"STRING API Error: {e.response.status_code} - {e.response.text[:100]}"
    except requests.exceptions.RequestException as e: return pd.DataFrame(), f"STRING API Request Failed: {e}"
    except pd.errors.EmptyDataError: return pd.DataFrame(), f"No data returned from STRING for {gene_symbol}."
    except ValueError as e: return pd.DataFrame(), f"Error parsing STRING data, columns might have changed: {e}"

@st.cache_data(ttl=3600*2)
def fetch_string_static_image_url(gene_symbol, organism, score_threshold, limit):
    species_id = ORGANISM_TO_TAXID.get(organism)
    if not species_id: return None, f"Organism '{organism}' not mapped for STRING static image."
    effective_limit = min(limit, 50)
    url = (f"https://string-db.org/api/image/network?identifiers={gene_symbol.upper()}"
           f"&species={species_id}&required_score={score_threshold}&limit={effective_limit}"
           f"&network_flavor=confidence&caller_identity=GenVedix_App_v1")
    return url, None

# --- Helper function to parse FASTA for GC content (specific to GC tool) ---
def parse_fasta_sequences_for_gc(fasta_string):
    sequences = []
    fasta_io = io.StringIO(fasta_string)
    for record in SeqIO.parse(fasta_io, "fasta"):
        seq_str = str(record.seq).upper()
        if record.seq and any(c in "ATGCU" for c in seq_str):
            sequences.append({"id": record.id, "seq": seq_str})
        elif record.seq:
            st.warning(f"Sequence '{record.id}' may not be DNA/RNA or contains many non-standard characters; skipped if unsuitable for GC analysis.")
    return sequences

# --- Helper function to calculate GC content profile ---
def calculate_gc_profile(sequence_str, window_size):
    positions = []
    gc_values = []
    sequence_str = sequence_str.upper().replace("U", "T")
    if not sequence_str or window_size <= 0 or window_size > len(sequence_str):
        return positions, gc_values
    for i in range(len(sequence_str) - window_size + 1):
        window = sequence_str[i:i + window_size]
        g_count = window.count('G')
        c_count = window.count('C')
        atgc_count = sum(window.count(b) for b in "ATGC")
        if atgc_count == 0 : gc_content = 0
        else: gc_content = (g_count + c_count) / float(atgc_count) if atgc_count >= window_size * 0.5 else 0

        positions.append(i + window_size // 2)
        gc_values.append(gc_content * 100)
    return positions, gc_values

# --- Helper function for nucleotide composition ---
def calculate_nucleotide_composition(dna_sequence):
    dna_sequence = dna_sequence.upper()
    counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0, 'U': 0, 'Other': 0}
    total_len = len(dna_sequence)
    if total_len == 0:
        return counts, total_len

    for base in dna_sequence:
        if base in counts:
            counts[base] += 1
        else:
            counts['Other'] += 1
    return counts, total_len

# --- Helper function for Protein Residue Counter ---
def calculate_amino_acid_composition(protein_sequence):
    protein_sequence = protein_sequence.strip().upper()

    standard_aas = "ACDEFGHIKLMNPQRSTVWY"
    final_counts = {aa: 0 for aa in standard_aas}
    final_counts["Other"] = 0 # For non-standard or ambiguous characters

    total_length = len(protein_sequence)
    if total_length == 0:
        data = [{"Amino Acid": aa, "Count": 0, "Percentage": 0.0} for aa in list(standard_aas) + ["Other"]]
        return pd.DataFrame(data)

    aa_counts_from_input = Counter(protein_sequence)

    for aa_char, count_val in aa_counts_from_input.items():
        if aa_char in final_counts:
            final_counts[aa_char] += count_val
        else:
            final_counts["Other"] += count_val

    composition_data = []
    for aa, count in final_counts.items():
        percentage = (count / total_length) * 100 if total_length > 0 else 0.0
        composition_data.append({"Amino Acid": aa, "Count": count, "Percentage": percentage})

    df = pd.DataFrame(composition_data)
    return df

# --- Helper function for Codon Usage ---
def calculate_codon_usage_from_cds(cds_sequences):
    if not cds_sequences:
        return pd.DataFrame(), {}, "No CDS sequences provided."

    standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
    codon_counts = Counter()
    aa_codon_counts = {aa: Counter() for aa in standard_table.protein_alphabet + "*"}
    total_codons = 0
    errors = []

    for i, seq_str in enumerate(cds_sequences):
        seq_str = seq_str.upper().replace("U", "T")
        if len(seq_str) % 3 != 0:
            errors.append(f"Sequence {i+1} (length {len(seq_str)}) is not a multiple of 3. Skipping.")
            continue

        if not all(c in "ATGC" for c in seq_str):
            non_dna_chars = set(c for c in seq_str if c not in "ATGC")
            errors.append(f"Sequence {i+1} contains non-DNA characters ({', '.join(sorted(list(non_dna_chars)))}). Skipping.")
            continue

        for j in range(0, len(seq_str), 3):
            codon = seq_str[j:j+3]
            codon_counts[codon] += 1
            total_codons += 1
            try:
                amino_acid = standard_table.forward_table.get(codon)
                if amino_acid:
                    aa_codon_counts[amino_acid][codon] += 1
                elif codon in standard_table.stop_codons:
                    aa_codon_counts["*"][codon] += 1
                else:
                    errors.append(f"Warning: Codon '{codon}' in sequence {i+1} is not a standard DNA codon.")
            except KeyError:
                errors.append(f"Warning: Codon '{codon}' in sequence {i+1} is not in the standard translation table.")

    if total_codons == 0:
        if not errors: return pd.DataFrame(), {}, "No valid codons found in the provided sequences to process."
        else: return pd.DataFrame(), {}, "No valid codons processed. Errors:\n" + "\n".join(errors)

    usage_data = []
    for codon, count in sorted(codon_counts.items()):
        frequency_permille = (count / total_codons) * 1000 if total_codons > 0 else 0
        amino_acid = standard_table.forward_table.get(codon)
        if amino_acid is None and codon in standard_table.stop_codons: amino_acid = "*"
        elif amino_acid is None: amino_acid = "?"
        usage_data.append({
            "Codon": codon, "Amino Acid": amino_acid,
            "Count": count, "Frequency (‰)": frequency_permille
        })

    df_usage = pd.DataFrame(usage_data)
    error_message_out = "\n".join(errors) if errors else None
    return df_usage, aa_codon_counts, error_message_out

# --- Helper functions for Phylogenetic Tree Viewer ---
@st.cache_data(max_entries=5, ttl=3600)
def parse_newick_tree(newick_string):
    try:
        handle = io.StringIO(newick_string)
        tree = Phylo.read(handle, "newick")
        return tree, None
    except Exception as e:
        return None, f"Error parsing Newick tree: {str(e)}. Please ensure valid Newick format."

def plot_phylo_tree_plotly(tree):
    if not tree or not tree.root: return None
    def get_clade_coords(clade, x_start=0, y_start=0, y_step=1, x_scale=1.0, coords=None, positions=None):
        if coords is None: coords = {'x': [], 'y': []}
        if positions is None: positions = {}
        current_x = x_start + (clade.branch_length * x_scale if clade.branch_length is not None else 0.5 * x_scale)
        if clade.is_terminal():
            positions[clade] = (current_x, y_start)
            coords['x'].extend([x_start, current_x, None]); coords['y'].extend([y_start, y_start, None])
            return y_start + y_step
        children = clade.clades
        if not children: return y_start
        child_y_positions = []; next_y = y_start
        for child in children:
            next_y = get_clade_coords(child, current_x, next_y, y_step, x_scale, coords, positions)
            child_y_positions.append(positions[child][1])
        current_node_y_mid = sum(child_y_positions) / len(child_y_positions) if child_y_positions else y_start
        positions[clade] = (current_x, current_node_y_mid)
        coords['x'].extend([x_start, current_x, None]); coords['y'].extend([current_node_y_mid, current_node_y_mid, None])
        if len(child_y_positions) > 1:
            min_child_y = min(child_y_positions); max_child_y = max(child_y_positions)
            coords['x'].extend([current_x, current_x, None]); coords['y'].extend([min_child_y, max_child_y, None])
        return next_y
    coordinates = {'x': [], 'y': []}; node_positions = {}
    num_terminals = len(tree.get_terminals()); y_step_val = 1.0
    get_clade_coords(tree.root, x_start=0, y_start=0, y_step=y_step_val, x_scale=1.0, coords=coordinates, positions=node_positions)
    node_x_coords = [pos[0] for pos in node_positions.values()]
    node_y_coords = [pos[1] for pos in node_positions.values()]
    node_text_labels = [clade.name if clade.name else (f"Internal_{idx}" if not clade.is_terminal() else "") for idx, clade in enumerate(node_positions.keys())]
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=coordinates['x'], y=coordinates['y'], mode='lines', line=dict(color='#4CAF50', width=1.5), hoverinfo='none'))
    fig.add_trace(go.Scatter(x=node_x_coords, y=node_y_coords, mode='markers+text', marker=dict(symbol='circle', size=6, color='#00bcd4'),
        text=node_text_labels, textposition="middle right", textfont=dict(color="#eaffea", size=10), hovertext=node_text_labels, hoverinfo='text'))
    fig.update_layout(title_text="Phylogenetic Tree", showlegend=False, xaxis=dict(title="Branch Length / Divergence", zeroline=False, showgrid=False),
        yaxis=dict(zeroline=False, showgrid=False, showticklabels=False, autorange="reversed"),
        plot_bgcolor="#23272f", paper_bgcolor="#23272f", font=dict(color="#eaffea"), height=max(400, num_terminals * 20 + 100))
    return fig

# --- Helper function for ORF Finder ---
def find_orfs(sequence_str, min_orf_len_bp=100, start_codons=["ATG"], stop_codons=["TAA", "TAG", "TGA"]):
    orfs = []
    seq_obj = Seq(sequence_str.upper())
    standard_table = CodonTable.unambiguous_dna_by_name["Standard"] # Standard NCBI translation table 1

    for strand_direction, nuc_sequence in [(1, seq_obj), (-1, seq_obj.reverse_complement())]:
        for frame_offset in range(3): # 0, 1, 2
            current_pos = frame_offset
            while current_pos + 2 < len(nuc_sequence):
                codon = str(nuc_sequence[current_pos : current_pos + 3])
                if codon in start_codons:
                    orf_start_in_frame = current_pos
                    orf_end_in_frame = -1
                    for j in range(orf_start_in_frame + 3, len(nuc_sequence) - 2, 3):
                        potential_stop_codon = str(nuc_sequence[j : j + 3])
                        if potential_stop_codon in stop_codons:
                            orf_end_in_frame = j + 2 # end of the stop codon
                            break
                    if orf_end_in_frame != -1:
                        orf_dna_seq_in_frame = nuc_sequence[orf_start_in_frame : orf_end_in_frame + 1]
                        orf_len_bp = len(orf_dna_seq_in_frame)

                        if orf_len_bp >= min_orf_len_bp:
                            protein_full = ""
                            try:
                                protein_full = str(orf_dna_seq_in_frame.translate(table=standard_table, cds=False, to_stop=False))
                            except CodonTable.TranslationError as e:
                                protein_full = f"[Translation Error: {e}]"

                            protein_seq = protein_full
                            if protein_full.endswith("*") and not protein_full.startswith("[Translation Error"):
                                protein_seq = protein_full[:-1]

                            if strand_direction == 1:
                                original_start = orf_start_in_frame
                                original_end = orf_end_in_frame
                                display_frame = frame_offset + 1
                            else:
                                original_start = len(seq_obj) - 1 - orf_end_in_frame
                                original_end = len(seq_obj) - 1 - orf_start_in_frame
                                display_frame = -(frame_offset + 1)

                            orfs.append({
                                "Frame": display_frame,
                                "Start (0-indexed)": original_start,
                                "End (0-indexed)": original_end,
                                "Length (bp)": orf_len_bp,
                                "DNA Sequence": str(orf_dna_seq_in_frame),
                                "Protein Sequence": protein_seq,
                                "Protein Length (aa)": len(protein_seq) if not protein_seq.startswith("[Translation Error") else 0
                            })
                current_pos += 3
    return sorted(orfs, key=lambda x: (x["Start (0-indexed)"], x["Frame"]))


# --- Function to design primers using Primer3 ---
def design_primers_with_primer3(seq_id, sequence, prod_min, prod_max, opt_size, min_size, max_size, opt_tm, min_tm, max_tm):
    if primer3 is None:
        return None, "Primer3 library is not installed."
    try:
        seq_args = {
            'SEQUENCE_ID': seq_id,
            'SEQUENCE_TEMPLATE': sequence,
        }
        global_args = {
            'PRIMER_OPT_SIZE': int(opt_size),
            'PRIMER_MIN_SIZE': int(min_size),
            'PRIMER_MAX_SIZE': int(max_size),
            'PRIMER_OPT_TM': float(opt_tm),
            'PRIMER_MIN_TM': float(min_tm),
            'PRIMER_MAX_TM': float(max_tm),
            'PRIMER_PRODUCT_SIZE_RANGE': [[int(prod_min), int(prod_max)]],
            'PRIMER_PICK_LEFT_PRIMER': 1,
            'PRIMER_PICK_INTERNAL_OLIGO': 0,
            'PRIMER_PICK_RIGHT_PRIMER': 1,
            'PRIMER_EXPLAIN_FLAG': 1,
            'PRIMER_NUM_RETURN': 5
        }
        results = primer3.bindings.designPrimers(seq_args, global_args)
        if results.get('PRIMER_PAIR_NUM_RETURNED', 0) == 0:
            error_message = "No primer pairs found. "
            explanations = []
            if 'PRIMER_LEFT_EXPLAIN' in results: explanations.append(f"Left: {results['PRIMER_LEFT_EXPLAIN']}")
            if 'PRIMER_RIGHT_EXPLAIN' in results: explanations.append(f"Right: {results['PRIMER_RIGHT_EXPLAIN']}")
            if 'PRIMER_PAIR_EXPLAIN' in results: explanations.append(f"Pair: {results['PRIMER_PAIR_EXPLAIN']}")
            if explanations: error_message += "Reasons: " + "; ".join(explanations)
            else: error_message += "Consider relaxing parameters or checking sequence quality."
            return results, error_message
        return results, None
    except Exception as e:
        if "found an N" in str(e).lower() or "too short" in str(e).lower() or "product size range" in str(e).lower():
             return None, f"Primer3 Input Error: {str(e)}. Please check sequence for N's, length, and product size parameters."
        return None, f"An error occurred during primer design: {str(e)}"


# --- Home Page Rendering ---
def render_home_page():
    st.markdown("<div style='height:3.5em;'></div>", unsafe_allow_html=True)
    st.markdown("<h1 style='color:#4CAF50; font-weight:900; text-shadow:none; text-align:center; font-size:5rem; letter-spacing:3px; margin-bottom:0.1em;'>GenVedix</h1>", unsafe_allow_html=True)
    st.markdown("<div style='text-align:center; margin-top:-1.2em; margin-bottom:0.2em;'><span style='color:#00bcd4; font-size:2.5rem; font-weight:700; letter-spacing:1.5px;'>An Aisle To Gene World</span></div>", unsafe_allow_html=True)
    st.markdown("<div style='text-align:center; margin-top:1.2em; margin-bottom:2em;'><span style='color:#eaffea; font-size:1.35rem; font-style:italic; font-family:Georgia, \"Times New Roman\", serif; font-weight:500;'>A modern, unified platform for gene annotation, networks,<br>pathways, and sequence analysis—built for researchers and students.</span></div>", unsafe_allow_html=True)
    st.markdown("---")

    st.markdown("<h2 class='home-section-header'>Explore GenVedix Features</h2>", unsafe_allow_html=True)
    st.markdown("""
    <p style='text-align:center; color: #e0e0e0; font-size:1.15em; margin-bottom: 2em; line-height:1.7;'>
    Welcome to GenVedix! Discover a comprehensive suite of bioinformatics tools designed to simplify your research and analysis.
    Below is a guide to our key features. Click on each to learn more about how it works and how to use it.
    </p>
    """, unsafe_allow_html=True)

    # --- Genomic Toolkit Features ---
    st.markdown("<h3 class='section-subheader' style='text-align:left; margin-left:0; font-size:1.7em !important;'> Genomic Toolkit</h3>", unsafe_allow_html=True)

    with st.container():
        st.markdown('<div class="home-tool-expander">', unsafe_allow_html=True)
        with st.expander(" Gene Annotation Search", expanded=False):
            st.markdown("""
            <p>Fetch detailed information about genes from the NCBI Entrez database.</p>
            <strong>How it works:</strong>
            <ul>
                <li>Enter a gene symbol (e.g., <code>TP53</code>, <code>BRCA1</code>).</li>
                <li>Select the organism of interest (e.g., <code>Homo sapiens</code>).</li>
                <li>Optionally, filter results by annotation status (All, Fully Annotated, Partially Annotated).</li>
                <li>Click "Search Gene Annotations".</li>
            </ul>
            <strong>Output:</strong>
            <ul>
                <li>Displays gene symbol, official name, aliases, function summary, chromosomal location, exon count, and NCBI Gene ID.</li>
                <li>Provides direct links to NCBI nucleotide IDs associated with the gene.</li>
                <li>Allows viewing of GenBank sequence details (preview, length, features) for associated nucleotide IDs.</li>
                <li>Enables downloading of the gene's FASTA sequence and annotation details in CSV, JSON, or TXT format.</li>
            </ul>
            <p><em>Navigate to "Genomic Toolkit" > " Gene Search" tab to use this feature.</em></p>
            """, unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)

    with st.container():
        st.markdown('<div class="home-tool-expander">', unsafe_allow_html=True)
        with st.expander(" Gene Interaction Network", expanded=False):
            st.markdown("""
            <p>Visualize gene interaction networks using data from the STRING DB, rendered with Pyvis for interactivity.</p>
            <strong>How it works:</strong>
            <ul>
                <li>Enter a gene symbol (e.g., <code>EGFR</code>).</li>
                <li>Select the target organism.</li>
                <li>Adjust advanced options like minimum interaction score (0-1000, higher is more stringent) and the approximate maximum number of interactors to display.</li>
                <li>Choose to enable/disable physics for network layout and show/hide edge scores.</li>
                <li>Click "Fetch & Visualize Network".</li>
            </ul>
            <strong>Output:</strong>
            <ul>
                <li>An interactive network graph where nodes are proteins/genes and edges represent interactions. The query gene is highlighted.</li>
                <li>A static network image directly from STRING DB is also displayed for a quick overview.</li>
                <li>Download options for the interactive HTML network and the static PNG image.</li>
                <li>A direct link to view the network on the STRING DB website.</li>
            </ul>
            <p><em>Navigate to "Genomic Toolkit" > "Gene Network" tab to use this feature.</em></p>
            """, unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)

    with st.container():
        st.markdown('<div class="home-tool-expander">', unsafe_allow_html=True)
        with st.expander("Pathway Analysis", expanded=False):
            st.markdown("""
            <p>Identify and explore biological pathways associated with a specific gene from KEGG and Reactome databases.</p>
            <strong>How it works:</strong>
            <ul>
                <li>Enter a gene symbol.</li>
                <li>Select the organism.</li>
                <li>Click "Fetch Pathways".</li>
            </ul>
            <strong>Output:</strong>
            <ul>
                <li>Lists pathways from KEGG and Reactome where the gene is involved, including Pathway ID, Name, and a direct URL to the database entry.</li>
                <li>For KEGG pathways, you can load and view the pathway diagram image directly within the app and download it.</li>
                <li>For Reactome pathways, a link to the SVG diagram is provided, and you can attempt to load/preview and download the SVG.</li>
            </ul>
            <p><em>Navigate to "Genomic Toolkit" > "Pathways" tab to use this feature.</em></p>
            """, unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)

    with st.container():
        st.markdown('<div class="home-tool-expander">', unsafe_allow_html=True)
        with st.expander("BLAST Search (NCBI)", expanded=False):
            st.markdown("""
            <p>Perform Basic Local Alignment Search Tool (BLAST) queries against NCBI databases.</p>
            <strong>How it works:</strong>
            <ul>
                <li>Paste your query sequence in FASTA format (or raw sequence, a default header will be added).</li>
                <li>Select the BLAST program (<code>blastn</code>, <code>blastp</code>, <code>blastx</code>, etc.).</li>
                <li>Choose the target NCBI database (e.g., <code>nr</code>, <code>nt</code>, <code>pdb</code>).</li>
                <li>Optionally, adjust the E-value threshold and maximum number of hits.</li>
                <li>Click "Run BLAST". (This may take some time depending on NCBI's servers and query complexity).</li>
            </ul>
            <strong>Output:</strong>
            <ul>
                <li>A table of significant hits, including Hit ID, Description, Score, E-value, Query Cover, Percent Identity, and Alignment Length.</li>
                <li>Download options for the BLAST results in CSV (tabular) and XML (raw NCBI output) formats.</li>
            </ul>
            <p><em>Navigate to "Genomic Toolkit" > "BLAST Search" tab to use this feature.</em></p>
            """, unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)

    with st.container():
        st.markdown('<div class="home-tool-expander">', unsafe_allow_html=True)
        with st.expander("Primer Design (Primer3)", expanded=False):
            st.markdown("""
            <p>Design PCR primers for a DNA sequence template using the Primer3 algorithm.</p>
            <strong>How it works:</strong>
            <ul>
                <li>Enter your DNA template sequence (raw or FASTA format).</li>
                <li>Specify parameters for primer design:
                    <ul>
                        <li>Product Size Range (min and max length of the PCR product).</li>
                        <li>Primer Size (optimal, min, max length of the primers).</li>
                        <li>Primer Melting Temperature (Tm) (optimal, min, max).</li>
                    </ul>
                </li>
                <li>Click "Design Primers".</li>
            </ul>
            <strong>Output:</strong>
            <ul>
                <li>A table listing suggested primer pairs, including their sequences, start positions, lengths, Tm, GC content, and expected product size.</li>
                <li>Option to download the primer results as a CSV file.</li>
                <li>An expander to view the full raw dictionary output from Primer3 for detailed diagnostics.</li>
            </ul>
            <p><em>Navigate to "Genomic Toolkit" > "Primer Design" tab to use this feature. Requires `primer3-py` and `primer3_core`.</em></p>
            """, unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)


    # --- FASTA Analysis Utilities Features ---
    st.markdown("<h3 class='section-subheader' style='text-align:left; margin-left:0; font-size:1.7em !important; margin-top:2em;'>FASTA Analysis Utilities</h3>", unsafe_allow_html=True)

    with st.container():
        st.markdown('<div class="home-tool-expander">', unsafe_allow_html=True)
        with st.expander("GC Profile Visualizer", expanded=False):
            st.markdown("""
            <p>Visualize the GC content distribution along DNA/RNA sequences.</p>
            <strong>How it works:</strong>
            <ul>
                <li>Paste one or more sequences in FASTA format, or upload a FASTA file.</li>
                <li>Set the "Window Size" (in base pairs/nucleotides) for calculating GC content in sliding windows.</li>
                <li>Click "Generate GC Profile Plot".</li>
            </ul>
            <strong>Output:</strong>
            <ul>
                <li>A line plot for each sequence (up to 5 are plotted for performance), showing GC percentage on the Y-axis and position (midpoint of the window) on the X-axis.</li>
                <li>Helps identify GC-rich or AT-rich regions within a sequence.</li>
            </ul>
            <p><em>Navigate to "FASTA Analysis Utilities" page to use this feature.</em></p>
            """, unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)

    with st.container():
        st.markdown('<div class="home-tool-expander">', unsafe_allow_html=True)
        with st.expander("Base Composition Analyzer", expanded=False):
            st.markdown("""
            <p>Calculate and visualize the nucleotide (A, T/U, G, C, Other) composition of a single DNA/RNA sequence.</p>
            <strong>How it works:</strong>
            <ul>
                <li>Paste a single raw DNA or RNA sequence (not FASTA).</li>
                <li>Click "Analyze Base Composition".</li>
            </ul>
            <strong>Output:</strong>
            <ul>
                <li>A pie chart showing the percentage of each nucleotide.</li>
                <li>A summary of the counts and percentages for each base present.</li>
            </ul>
            <p><em>Navigate to "FASTA Analysis Utilities" page to use this feature.</em></p>
            """, unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)

    with st.container():
        st.markdown('<div class="home-tool-expander">', unsafe_allow_html=True)
        with st.expander("ORF Finder", expanded=False):
            st.markdown("""
            <p>Identify Open Reading Frames (ORFs) within a DNA sequence across all six reading frames.</p>
            <strong>How it works:</strong>
            <ul>
                <li>Upload a FASTA file (the first sequence will be used) or paste a raw DNA sequence.</li>
                <li>Specify the "Minimum ORF length" in amino acids (default is 30 aa).</li>
                <li>The tool uses standard start (ATG) and stop (TAA, TAG, TGA) codons.</li>
                <li>Click "Find ORFs".</li>
            </ul>
            <strong>Output:</strong>
            <ul>
                <li>A table listing each found ORF, including its frame, start/end positions (0-indexed), length in base pairs, translated protein length, and a preview of the protein sequence.</li>
                <li>Download buttons for the DNA sequences of all found ORFs and their corresponding protein sequences, both in FASTA format.</li>
            </ul>
            <p><em>Navigate to "FASTA Analysis Utilities" page to use this feature.</em></p>
            """, unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)

    with st.container():
        st.markdown('<div class="home-tool-expander">', unsafe_allow_html=True)
        with st.expander("3D Structure Viewer (PDB)", expanded=False):
            st.markdown("""
            <p>Visualize 3D macromolecular structures directly by fetching data from the Protein Data Bank (PDB).</p>
            <strong>How it works:</strong>
            <ul>
                <li>Enter a 4-character PDB ID (e.g., <code>1MOL</code>, <code>6LU7</code>).</li>
                <li>Adjust the viewer height and select preferred display style (e.g., cartoon, stick, sphere) and color scheme.</li>
                <li>Click "Visualize Structure".</li>
            </ul>
            <strong>Output:</strong>
            <ul>
                <li>An interactive 3D view of the protein structure, rendered using py3Dmol. You can rotate, zoom, and pan the structure with your mouse.</li>
            </ul>
            <p><em>Navigate to "FASTA Analysis Utilities" page to use this feature. Requires `py3Dmol` library.</em></p>
            """, unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)

    with st.container():
        st.markdown('<div class="home-tool-expander">', unsafe_allow_html=True)
        with st.expander("Protein Residue Counter", expanded=False):
            st.markdown("""
            <p>Calculate and display the frequency of each amino acid in a protein sequence.</p>
            <strong>How it works:</strong>
            <ul>
                <li>Paste a single raw protein sequence (not FASTA).</li>
                <li>Click "Calculate Amino Acid Frequencies".</li>
            </ul>
            <strong>Output:</strong>
            <ul>
                <li>A bar chart showing the percentage frequency of each amino acid.</li>
                <li>A table detailing the count and percentage for each of the 20 standard amino acids and any 'Other' non-standard residues.</li>
            </ul>
            <p><em>Navigate to "FASTA Analysis Utilities" page to use this feature.</em></p>
            """, unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)

    with st.container():
        st.markdown('<div class="home-tool-expander">', unsafe_allow_html=True)
        with st.expander("Expression Pattern Mapper (Heatmap)", expanded=False):
            st.markdown("""
            <p>Generate a heatmap to visualize expression patterns from tabular data (e.g., gene expression across samples).</p>
            <strong>How it works:</strong>
            <ul>
                <li>Upload a CSV file where rows typically represent genes/features and columns represent samples/conditions.</li>
                <li>Select the column to use as the index (genes/features).</li>
                <li>Select the data columns to include in the heatmap (samples/conditions). These must be numeric.</li>
                <li>Choose a colorscale for the heatmap.</li>
                <li>Optionally, apply Z-score normalization either by row (gene-wise) or by column (sample-wise).</li>
                <li>Click "Generate Expression Heatmap".</li>
            </ul>
            <strong>Output:</strong>
            <ul>
                <li>An interactive heatmap displaying the data, with options to hover for values.</li>
                <li>An expander to view the processed data table that was used to generate the heatmap.</li>
            </ul>
            <p><em>Navigate to "FASTA Analysis Utilities" page to use this feature.</em></p>
            """, unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)

    with st.container():
        st.markdown('<div class="home-tool-expander">', unsafe_allow_html=True)
        with st.expander("Codon Usage Analyzer", expanded=False):
            st.markdown("""
            <p>Analyze codon usage patterns from one or more coding DNA sequences (CDS).</p>
            <strong>How it works:</strong>
            <ul>
                <li>Upload a FASTA file containing CDS sequences or paste raw CDS sequences (one per line, no headers).</li>
                <li>The tool assumes standard DNA codons.</li>
                <li>Click "Calculate Codon Usage".</li>
            </ul>
            <strong>Output:</strong>
            <ul>
                <li>A table showing each codon, the amino acid it codes for, its count across all input sequences, and its frequency per mille (‰).</li>
                <li>Bar charts for each amino acid (excluding '?'), displaying the counts of its synonymous codons, helping to visualize codon preference.</li>
            </ul>
            <p><em>Navigate to "FASTA Analysis Utilities" page to use this feature.</em></p>
            """, unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)

    with st.container():
        st.markdown('<div class="home-tool-expander">', unsafe_allow_html=True)
        with st.expander("Phylogenetic Tree Viewer (Newick)", expanded=False):
            st.markdown("""
            <p>Visualize phylogenetic trees provided in the Newick file format.</p>
            <strong>How it works:</strong>
            <ul>
                <li>Upload a tree file with a <code>.nwk</code>, <code>.newick</code>, <code>.tree</code>, or <code>.txt</code> extension containing a valid Newick string.</li>
            </ul>
            <strong>Output:</strong>
            <ul>
                <li>An interactive plot of the phylogenetic tree, rendered using Plotly. You can see branch lengths and node labels.</li>
                <li>An expander to view the raw Newick string that was parsed.</li>
            </ul>
            <p><em>Navigate to "FASTA Analysis Utilities" page to use this feature.</em></p>
            """, unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)


    # --- Metagenomics Resources ---
    st.markdown("<h3 class='section-subheader' style='text-align:left; margin-left:0; font-size:1.7em !important; margin-top:2em;'>Metagenomics Resources</h3>", unsafe_allow_html=True)
    with st.container():
        st.markdown('<div class="home-tool-expander">', unsafe_allow_html=True)
        with st.expander("External Metagenomics Databases & Tools", expanded=False):
            st.markdown("""
            <p>Discover a curated collection of external web servers, databases, and command-line tools essential for metagenomic research.</p>
            <strong>How it works:</strong>
            <ul>
                <li>Navigate to the "Metagenomics Resources" page using the sidebar.</li>
                <li>Browse through categorized lists of tools, including general analysis platforms (e.g., MG-RAST, MGnify), taxonomic profilers (e.g., MetaPhlAn, Kraken2), functional annotators (e.g., eggNOG-mapper), assemblers, and specialized databases (e.g., CARD for antibiotic resistance).</li>
                <li>Each tool has a brief description and a button to "Open Tool".</li>
                <li>For tools marked as "iframe-able", clicking the button will attempt to embed the external tool directly within GenVedix for convenience. Otherwise, it will open in a new browser tab.</li>
            </ul>
            <p><em>This section serves as a quick-access directory to popular resources in the metagenomics field.</em></p>
            """, unsafe_allow_html=True)
        st.markdown('</div>', unsafe_allow_html=True)

    st.markdown("---")
    st.markdown("<h3 style='color: #00bcd4; text-align:center; margin-top:2em; margin-bottom:1em;'>Getting Started with GenVedix</h3>", unsafe_allow_html=True)
    st.markdown("""
    <p style='text-align:center; color: #eaffea; font-size:1.1em;'>
    To begin, simply use the sidebar navigation on the left to select a page and explore the tools:
    <ul style='text-align:left; max-width:550px; margin:auto; list-style-position: inside;'>
        <li><b>Home:</b> You are here! Get an overview of GenVedix and its features.</li>
        <li><b>Genomic Toolkit:</b> Your primary workbench for gene-centric analyses including search, networks, pathways, BLAST, and primer design.</li>
        <li><b>FASTA Analysis Utilities:</b> A collection of tools for sequence manipulation, visualization, and property calculation.</li>
        <li><b>Metagenomics Resources:</b> A directory of valuable external tools and databases for metagenomic studies.</li>
        <li><b>About:</b> Learn more about the project, the developer, and acknowledgements.</li>
    </ul>
    </p>
    """, unsafe_allow_html=True)
    st.markdown("<div style='height:2em;'></div>", unsafe_allow_html=True) # Add some bottom spacing


# ---About Section ---
def render_about_section():
    st.markdown("<h1 class='main-header'>About GenVedix</h1>", unsafe_allow_html=True)
    with st.container():
        st.markdown("""
        <div class='main-container'>
            <h2 style='color: #4CAF50!important; border-bottom: 2px solid #4CAF50; padding-bottom: 0.5rem;'>Developer's Desk</h2>
            <div style='display: flex; flex-direction: column; align-items: center; margin: 2rem 0;'>
                <img src='https://media.licdn.com/dms/image/v2/D5603AQFfI1KWVSWl0Q/profile-displayphoto-shrink_400_400/B56ZRRpRfTH0Ag-/0/1736536563970?e=1752710400&v=beta&t=4UtS-OUEjdd0dV7BMRQFA2CxjhOd-tgFPmYbeNAzqB4'
                     style='width: 180px; height: 180px; object-fit: cover; border-radius: 50%; border: 3px solid #4CAF50; margin-bottom: 1rem;'>
                <h3>Tejal Kale</h3>
                <p style='color: #4CAF50!important; margin-bottom: 0.5rem;'>Developer of GenVedix</p>
                <div class='feature-card' style='max-width: 800px;'>
                    <p style='text-align: justify; line-height: 1.6;'>
                    Hello Bioinformaticians! I’m Tejal Kale, a Master’s student in Bioinformatics at Deccan Education Society’s Pune University.

GenVedix is more than just a mini-project, it’s a solution born from lived experience. During my journey through gene exploration, I often found myself overwhelmed by fragmented tools, scattered databases, and time-consuming processes. Like many of you, I longed for a unified, accessible platform one that could simplify my workflow and help me focus on understanding genes, not chasing them across multiple sources.

That’s where GenVedix was born. It is built for you, the curious student trying to decode your first gene, the researcher seeking quick and reliable insights, the educator wanting to demonstrate genome structure clearly, and anyone passionate about unraveling biological complexity. GenVedix combines essential gene analysis functions into one cohesive interface clean, intuitive, and beginner friendly, with no installation needed.

I developed GenVedix not just as a tool, but as a companion in your genomic journey. Whether you're analyzing gene features, comparing sequences, or simply learning the ropes, GenVedix is designed to feel like an extension of your curiosity reliable, insightful, and always ready to assist.
                    </p>
                </div>
            </div>
        </div>
        """, unsafe_allow_html=True)
    with st.container():
        st.markdown("""
        <div class='main-container' style='margin-top: 2rem;'>
            <h2 style='color: #4CAF50!important; border-bottom: 2px solid #4CAF50; padding-bottom: 0.5rem;'>Meet My Mentor</h2>
            <div style='display: flex; flex-direction: column; align-items: center; margin: 2rem 0;'>
                <img src='https://despu.edu.in/assets/admin/img/faculty_images/1738156003_375696d99642e3a89742.jpg'
                     style='width: 180px; height: 180px; object-fit: cover; border-radius: 50%; border: 3px solid #4CAF50; margin-bottom: 1rem;'>
                <h3>Dr. Kushagra Kashyap</h3>
                <p style='color: #4CAF50!important; margin-bottom: 0.5rem;'>Assistant Professor - DES Pune University</p>
                <div class='feature-card' style='max-width: 800px;'>
                    <p style='text-align: justify; line-height: 1.6;'>
                    Meet Dr. Kushagra Kashyap, a passionate and enthusiastic educator with a PhD from CSIR-CDRI Lucknow.
                    As an Assistant Professor at Deccan Education Society's Pune University, Dr. Kashyap brings his expertise
                    in Structural Bioinformatics and Cheminformatics to the classroom, inspiring students to explore the
                    fascinating world of computational biology and drug discovery. With a strong commitment to academic excellence,
                    Dr. Kashyap's teaching style is characterized by clarity, enthusiasm, and a genuine passion for mentoring
                    the next generation of scientists.
                    </p>
                </div>
            </div>
        </div>
        """, unsafe_allow_html=True)
    with st.container():
        st.markdown("""
        <div class='main-container' style='margin-top: 2rem;'>
            <h2 style='color: #4CAF50!important; border-bottom: 2px solid #4CAF50; padding-bottom: 0.5rem;'>Acknowledgements</h2>
            <div class='feature-card' style='margin: 2rem 0; padding: 1.5rem;'>
                <p style='color: #e0e0e0; line-height: 1.6; font-size: 16px;'>
                I would like to extend my heartfelt gratitude to:<br><br>
                1. Dr. Kushagra Kashyap, for his invaluable guidance and support throughout this project. His expertise and enthusiasm were instrumental in shaping my work.<br><br>
                2. Dr. Poonam Deshpande, for her exceptional guidance, mentorship, and unwavering support. Her insights and feedback were crucial to the project's success.<br><br>
                3. DES Pune University, for providing me with the necessary resources, infrastructure, and opportunities to complete this project.
                </p>
            </div>
        </div>
        """, unsafe_allow_html=True)
    with st.container():
        st.markdown("<div class='main-container' style='margin-top: 2rem;'>", unsafe_allow_html=True)
        st.markdown("<h2 style='color: #4CAF50!important; border-bottom: 2px solid #4CAF50; padding-bottom: 0.5rem;'>Purpose & Overview</h2>", unsafe_allow_html=True)
        st.markdown("""
        <div class='title-box' style='background-color: #2a2a2a; padding: 1.5rem; border-radius: 8px; margin-bottom: 1.5rem; border-left: 5px solid #4CAF50;'>
            <h3 style='color: #4CAF50!important; margin-bottom: 0.5rem;'>GenVedix: Your Comprehensive Gene Exploration Platform</h3>
            <p style='color: #e0e0e0;'>GenVedix is designed to empower students, educators, and early-stage researchers by offering a streamlined, user-friendly platform that brings together essential genomic tools in one place. From real-time gene annotation and sequence analysis to GC content calculation, codon usage insights, ORF detection, and pairwise alignment, GenVedix simplifies complex bioinformatics tasks. Integrated with trusted databases like NCBI and built on a clean, accessible interface, it enables users to explore, interpret, and visualize gene-related data without technical barriers—making learning and research more intuitive, efficient, and engaging.</p>
        </div>
        """, unsafe_allow_html=True)
        st.markdown("""
            <div class='feature-card' style='margin-bottom: 1.5rem; padding: 1.5rem;'>
                <h4 style='color: #4CAF50!important;'>Core Objectives</h4>
                <p style='color: #e0e0e0; line-height: 1.6;'>
                GenVedix was developed to streamline bioinformatic workflows, providing researchers with:
                <ul style='color: #e0e0e0; padding-left: 20px; margin-top: 0.5rem;'>
                    <li>A comprehensive <strong>Genomic Toolkit</strong> including Gene Search, Network Analysis, Pathway Exploration, BLAST, and Primer Design.</li>
                    <li>A versatile suite of <strong>FASTA Analysis Utilities</strong> for GC Content Profiling, Base Composition, ORF Finding, 3D Structure Visualization (PDB), Protein Residue Counting, Expression Heatmapping, Codon Usage Analysis, and Phylogenetic Tree Viewing.</li>
                    <li>A curated list of external <strong>Metagenomics Resources</strong>.</li>
                    <li>A user-friendly interface designed for both experts and students.</li>
                    <li>Centralized access to data from trusted sources like NCBI Entrez, STRING DB, KEGG, Reactome, and RCSB PDB.</li>
                </ul>
                The tool aims to reduce time spent on data gathering and tool switching, thereby increasing focus on analysis and discovery.
                </p>
            </div>
        """, unsafe_allow_html=True)
        st.markdown("</div>", unsafe_allow_html=True)
    with st.container():
        st.markdown("""
        <div class='main-container' style='margin-top: 2rem;'>
            <h2 style='color: #4CAF50!important; border-bottom: 2px solid #4CAF50; padding-bottom: 0.5rem;'>Contact & Feedback</h2>
<div class='feature-card' style='margin: 2rem 0; padding: 1.5rem;'>
    <h3 style='color: #4CAF50!important;'>Tejal Kale</h3>
    <p style='color: #e0e0e0 ; line-height: 1.6; font-size: 16px;'>
        Master's Student in Bioinformatics<br>
        DES Pune University<br>
        <span style='color: #4CAF50;'>📧 </span><a href="mailto:kaletejal05@mail.com" style='color: #4CAF50;'>kaletejal05@mail.com</a>
    </p>
    <div style='margin-top: 1rem;'>
        <a href='https://www.linkedin.com/in/tejal-kale-0b2195299?utm_source=share&utm_campaign=share_via&utm_content=profile&utm_medium=android_app' target='_blank' style='text-decoration: none; margin-right: 1rem;'>
            <button class='custom-link-button'>
                Connect on LinkedIn
            </button>
        </a>
        <a href='https://github.com/TejalKale05' target='_blank' style='text-decoration: none;'>
            <button class='custom-link-button'>
                View GitHub Profile
            </button>
        </a>
    </div>
</div>
<div class='feature-card' style='margin: 2rem 0; padding: 1.5rem;'>
    <h3 style='color: #4CAF50!important;'>Dr. Kushagra Kashyap</h3>
    <p style='color: #e0e0e0; line-height: 1.6; font-size: 16px;'>
        Assistant Professor<br>
        DES Pune University
    </p>
    <div style='margin-top: 1rem;'>
        <a href='https://www.linkedin.com/in/dr-kushagra-kashyap-b230a3bb' target='_blank' style='text-decoration: none;'>
            <button class='custom-link-button'>
                Connect on LinkedIn
            </button>
        </a>
    </div>
</div>
""", unsafe_allow_html=True)
        st.markdown("</div>", unsafe_allow_html=True)

# --- Sections for Genomic Toolkit Page ---
def render_gene_search_section():
    st.subheader("Gene Annotation Search")
    st.markdown("<p style='font-size: 0.9em; color: #b0b0b0;'>Search for detailed gene annotations using NCBI Entrez.</p>", unsafe_allow_html=True)
    with st.expander("Gene Search Filters", expanded=True):
        database_option_gs = st.radio("Annotation Status Filter", ["All", "Fully Annotated Only", "Partially Annotated"], key="db_filter_gene_search_section_gt")
        selected_organism_gs = st.selectbox("Select Organism", ORGANISM_OPTIONS, key="organism_gene_search_section_gt")
    gene_name_input_gs = st.text_input(" Enter gene name (e.g., INS, BRCA1):", key="gene_input_gs_gt")
    for k in ['gs_gene_annotations_gt', 'gs_current_gene_search_term_gt', 'gs_current_organism_search_gt', 'gs_gene_error_gt']:
        if k not in st.session_state: st.session_state[k] = None
    if st.button("Search Gene Annotations", key="search_button_gs_gt", use_container_width=True, type="primary"):
        st.session_state.gs_gene_annotations_gt = None; st.session_state.gs_gene_error_gt = None
        if gene_name_input_gs:
            with st.spinner(f"Searching for '{gene_name_input_gs}'..."):
                for key_s in list(st.session_state.keys()):
                    if key_s.startswith(("gs_seq_details_gt_", "gs_fasta_data_gt_", "gs_dl_data_gt_", "gs_dl_fname_gt_", "gs_dl_mime_gt_")):
                        del st.session_state[key_s]
                annotations, error = fetch_gene_annotation(gene_name_input_gs, selected_organism_gs)
                st.session_state.gs_gene_annotations_gt = annotations; st.session_state.gs_gene_error_gt = error
                st.session_state.gs_current_gene_search_term_gt = gene_name_input_gs; st.session_state.gs_current_organism_search_gt = selected_organism_gs
        else: st.warning("Please enter a gene name."); st.session_state.gs_gene_annotations_gt = None
    if st.session_state.gs_gene_error_gt: st.error(st.session_state.gs_gene_error_gt)
    if st.session_state.gs_gene_annotations_gt is not None and not st.session_state.gs_gene_error_gt:
        gene_annotations_to_display = st.session_state.gs_gene_annotations_gt
        if database_option_gs == "Fully Annotated Only": gene_annotations_to_display = [g for g in gene_annotations_to_display if g["Annotation Status"] == "Complete"]
        elif database_option_gs == "Partially Annotated": gene_annotations_to_display = [g for g in gene_annotations_to_display if g["Annotation Status"] == "Partial"]
        if not gene_annotations_to_display:
            if st.session_state.gs_current_gene_search_term_gt: st.warning(f"⚠️ No genes matching '{database_option_gs}' criteria found for '{st.session_state.gs_current_gene_search_term_gt}'.")
        else:
            st.success(f"✅ Found {len(gene_annotations_to_display)} result(s) for '{st.session_state.gs_current_gene_search_term_gt}' (Filter: {database_option_gs}).")
            def display_gene_details_actions_gs(gene_data, key_prefix="gs_gt_"):
                st.markdown("<div class='result-box'>", unsafe_allow_html=True)
                col1_info, col2_info = st.columns([3,1])
                with col1_info:
                    st.markdown(f"<h6>{gene_data.get('Gene Symbol', 'N/A')} - {gene_data.get('Official Name', 'N/A')}</h6>", unsafe_allow_html=True)
                    st.markdown(f"**Aliases:** {gene_data.get('Aliases', 'N/A')}"); st.markdown(f"**Location:** Chromosome {gene_data.get('Chromosome', 'N/A')}, Position {gene_data.get('Position', 'N/A')}")
                    st.markdown(f"**Exons:** {gene_data.get('Exons', 'N/A')}"); st.markdown(f"**Function:** <div style='max-height: 150px; overflow-y: auto; font-size: 0.9em;'>{gene_data.get('Function', 'N/A')}</div>", unsafe_allow_html=True)
                with col2_info:
                    st.markdown(f"**Organism:** {gene_data.get('Source Organism', 'N/A')}"); st.markdown(f"**Gene ID:** {gene_data.get('Gene ID', 'N/A')}")
                    st.markdown(f"**Annotation:** {gene_data.get('Annotation Status', 'N/A')}"); st.markdown(f"**NCBI Nucleotide IDs:** {gene_data.get('Nucleotide IDs', 'N/A')}")
                st.markdown("---")
                action_col1, action_col2 = st.columns(2)
                nuccore_id_to_fetch = gene_data.get('Nucleotide IDs', "N/A").split(", ")[0] if gene_data.get('Nucleotide IDs', "N/A") != "N/A" else None
                with action_col1:
                    st.markdown("<h6>Sequence Details</h6>", unsafe_allow_html=True)
                    if nuccore_id_to_fetch:
                        seq_details_key = f"{key_prefix}seq_details_{gene_data.get('Gene ID', 'unknown')}"
                        if st.button("View Sequence", key=f"{key_prefix}view_seq_{gene_data.get('Gene ID', 'unknown')}", help="Fetch/display GenBank details"):
                            with st.spinner(f"Fetching sequence for {nuccore_id_to_fetch}..."): seq_data, seq_err = fetch_gene_sequence(nuccore_id_to_fetch); st.session_state[seq_details_key] = (seq_data, seq_err)
                        if seq_details_key in st.session_state:
                            seq_data_gs, seq_err_gs = st.session_state[seq_details_key]
                            if seq_err_gs: st.error(f"Sequence Error: {seq_err_gs}")
                            if seq_data_gs:
                                with st.expander("GenBank Details", expanded=False):
                                    st.markdown(f"**NCBI ID:** {seq_data_gs.get('id', 'N/A')}"); st.markdown(f"**Description:** {seq_data_gs.get('description', 'N/A')}")
                                    st.markdown(f"**Length:** {seq_data_gs.get('length', 'N/A')} bp"); st.text_area("Sequence Preview (first 1000bp)", value=seq_data_gs.get('sequence', ''), height=150, disabled=True, key=f"{key_prefix}text_seq_preview_{gene_data.get('Gene ID', 'unknown')}")
                                    if seq_data_gs.get('features'):
                                        st.markdown("<span style='color:#8bc34a; font-weight:700;'>Features (first 10):</span>", unsafe_allow_html=True)
                                        for feat in seq_data_gs['features'][:10]: st.markdown(f"<span style='color:#eaffea;'>• <b>Type:</b> {feat.get('type', 'N/A')}, <b>Location:</b> {feat.get('location', 'N/A')}</span>", unsafe_allow_html=True)
                    else: st.caption("No Nucleotide ID for sequence view.")
                with action_col2:
                    st.markdown("<h6>FASTA Sequence</h6>", unsafe_allow_html=True)
                    if nuccore_id_to_fetch:
                        fasta_data_key = f"{key_prefix}fasta_data_{gene_data.get('Gene ID', 'unknown')}"
                        if st.button("Prepare FASTA", key=f"{key_prefix}dl_fasta_btn_{gene_data.get('Gene ID', 'unknown')}", help="Prepare FASTA for download"):
                            with st.spinner(f"Generating FASTA for {nuccore_id_to_fetch}..."): fasta_content, fasta_err = fetch_fasta_sequence_data(nuccore_id_to_fetch); st.session_state[fasta_data_key] = (fasta_content, fasta_err)
                        if fasta_data_key in st.session_state:
                            fasta_content_gs, fasta_err_gs = st.session_state[fasta_data_key]
                            if fasta_err_gs: st.error(f"FASTA Error: {fasta_err_gs}")
                            if fasta_content_gs: st.download_button("⬇️ Download FASTA", fasta_content_gs, f"{gene_data.get('Gene Symbol', 'gene')}_{nuccore_id_to_fetch}.fasta", "text/plain", key=f"{key_prefix}dl_fasta_final_{gene_data.get('Gene ID', 'unknown')}")
                    else: st.caption("No Nucleotide ID for FASTA.")
                st.markdown("---"); st.markdown("<h6>Download Annotation Data</h6>", unsafe_allow_html=True)
                dl_fields_gs = st.multiselect("Select fields for download:", options=list(gene_data.keys()), default=["Gene Symbol", "Official Name", "Function", "Chromosome", "Position", "Exons", "Source Organism", "Gene ID"], key=f"{key_prefix}fields_dl_{gene_data.get('Gene ID', 'unknown')}")
                dl_format_gs = st.radio("Select download format:", ["CSV", "JSON", "TXT"], horizontal=True, key=f"{key_prefix}format_dl_{gene_data.get('Gene ID', 'unknown')}")
                dl_data_key = f"{key_prefix}dl_data_{gene_data.get('Gene ID', 'unknown')}"; dl_fname_key = f"{key_prefix}dl_fname_{gene_data.get('Gene ID', 'unknown')}"; dl_mime_key = f"{key_prefix}dl_mime_{gene_data.get('Gene ID', 'unknown')}"
                if st.button("Prepare Download", key=f"{key_prefix}prep_dl_{gene_data.get('Gene ID', 'unknown')}", type="secondary"):
                    filtered_data_dl_gs = {k: gene_data[k] for k in dl_fields_gs if k in gene_data}
                    if dl_format_gs == "CSV": st.session_state[dl_data_key] = pd.DataFrame([filtered_data_dl_gs]).to_csv(index=False); st.session_state[dl_fname_key] = f"{gene_data.get('Gene Symbol', 'gene')}_annotation.csv"; st.session_state[dl_mime_key] = "text/csv"
                    elif dl_format_gs == "JSON": st.session_state[dl_data_key] = json.dumps(filtered_data_dl_gs, indent=2); st.session_state[dl_fname_key] = f"{gene_data.get('Gene Symbol', 'gene')}_annotation.json"; st.session_state[dl_mime_key] = "application/json"
                    else: st.session_state[dl_data_key] = "\n".join([f"{k}: {v}" for k, v in filtered_data_dl_gs.items()]); st.session_state[dl_fname_key] = f"{gene_data.get('Gene Symbol', 'gene')}_annotation.txt"; st.session_state[dl_mime_key] = "text/plain"
                    st.success(f"{dl_format_gs} file prepared.")
                if dl_data_key in st.session_state: st.download_button(f"⬇️ Download {st.session_state[dl_fname_key]}", st.session_state[dl_data_key], st.session_state[dl_fname_key], st.session_state[dl_mime_key], key=f"{key_prefix}actual_dl_btn_{gene_data.get('Gene ID', 'unknown')}")
                st.markdown("</div>", unsafe_allow_html=True)
            if len(gene_annotations_to_display) > 1:
                tab_titles_gs = [f"{g.get('Gene Symbol', 'N/A')} ({g.get('Gene ID', 'N/A')})" for g in gene_annotations_to_display]; tabs_gs = st.tabs(tab_titles_gs)
                for i, tab_content_gs in enumerate(tabs_gs):
                    with tab_content_gs: display_gene_details_actions_gs(gene_annotations_to_display[i], key_prefix=f"gs_tab_gt_{i}_")
            elif len(gene_annotations_to_display) == 1: display_gene_details_actions_gs(gene_annotations_to_display[0], key_prefix="gs_single_gt_")
    elif st.session_state.gs_current_gene_search_term_gt and not st.session_state.gs_gene_error_gt: st.info(f"No results for '{st.session_state.gs_current_gene_search_term_gt}'.")

def render_gene_network_section():
    st.subheader("Gene Interaction Network (STRING DB with Pyvis)")
    st.markdown("<p style='font-size: 0.9em; color: #b0b0b0;'>Visualize interactive gene networks.</p>", unsafe_allow_html=True)
    col1, col2 = st.columns([2,1]);
    with col1: gene_symbol_gn = st.text_input("Enter Gene Symbol (e.g., TP53):", key="gene_symbol_gn_pyvis_gt").upper()
    with col2: selected_organism_gn = st.selectbox("Select Organism:", ORGANISM_OPTIONS, key="organism_gn_pyvis_gt", index=0)
    with st.expander("Advanced Network Options"):
        score_threshold_gn = st.slider("Min. Interaction Score:", 0, 1000, 400, 50, key="score_thresh_gn_pyvis_gt")
        limit_interactions_gn = st.slider("Max. Interactors/Edges (approx.):", 5, 50, 15, 1, key="limit_interact_gn_pyvis_gt")
        physics_gn = st.checkbox("Enable Physics Engine", value=True, key="physics_gn_pyvis_gt")
        show_edge_labels_gn = st.checkbox("Show Edge Scores", value=False, key="edge_labels_gn_pyvis_gt")
    for k in ['pyvis_html_gn_gt', 'static_string_image_url_gn_gt', 'string_image_bytes_gn_gt', 'string_image_filename_gn_gt', 'interactions_df_gn_gt', 'error_gn_gt']:
        if k not in st.session_state: st.session_state[k] = None
    if st.button("Fetch & Visualize Network", key="fetch_viz_network_gn_pyvis_gt", use_container_width=True, type="primary"):
        st.session_state.pyvis_html_gn_gt = None; st.session_state.static_string_image_url_gn_gt = None; st.session_state.string_image_bytes_gn_gt = None
        st.session_state.interactions_df_gn_gt = None; st.session_state.error_gn_gt = None
        if gene_symbol_gn and selected_organism_gn:
            with st.spinner(f"Fetching interaction data for {gene_symbol_gn}..."):
                interactions_df, error_interact = fetch_string_interaction_data_for_pyvis(gene_symbol_gn, selected_organism_gn, int(score_threshold_gn), int(limit_interactions_gn))
                st.session_state.interactions_df_gn_gt = interactions_df
                if error_interact: st.session_state.error_gn_gt = error_interact
            with st.spinner("Fetching static image URL..."):
                static_url, error_static_url = fetch_string_static_image_url(gene_symbol_gn, selected_organism_gn, int(score_threshold_gn), int(limit_interactions_gn))
                st.session_state.static_string_image_url_gn_gt = static_url
                if error_static_url and not st.session_state.error_gn_gt: st.session_state.error_gn_gt = error_static_url
            if st.session_state.error_gn_gt: st.error(st.session_state.error_gn_gt)
            if st.session_state.interactions_df_gn_gt is not None and not st.session_state.interactions_df_gn_gt.empty:
                with st.spinner("Generating interactive network..."):
                    G = nx.Graph()
                    for _, row in st.session_state.interactions_df_gn_gt.iterrows():
                        p1, p2, score = str(row['protein1']).upper(), str(row['protein2']).upper(), row['score']
                        G.add_node(p1, title=p1, size=10, color='#69b3a2'); G.add_node(p2, title=p2, size=10, color='#69b3a2')
                        G.add_edge(p1, p2, value=score/1000.0, title=f"Score: {score}" if show_edge_labels_gn else str(score), weight=score)
                    query_gene_upper = gene_symbol_gn.upper()
                    if query_gene_upper in G: G.nodes[query_gene_upper]['color'] = '#FF4B4B'; G.nodes[query_gene_upper]['size'] = 20
                    elif query_gene_upper: G.add_node(query_gene_upper, title=query_gene_upper, size=20, color='#FF4B4B')
                    net = Network(notebook=True, height="750px", width="100%", cdn_resources='remote', heading=f"Network for {gene_symbol_gn}"); net.from_nx(G)
                    if physics_gn: net.force_atlas_2based(gravity=-80, central_gravity=0.01, spring_length=100, spring_strength=0.1, damping=0.5)
                    else: net.options.physics.enabled = False
                    temp_html_file_path = f"temp_network_gt_{gene_symbol_gn}_{int(time.time())}.html"
                    try:
                        net.save_graph(temp_html_file_path)
                        if os.path.exists(temp_html_file_path):
                            with open(temp_html_file_path, "r", encoding="utf-8") as f: st.session_state.pyvis_html_gn_gt = f.read()
                            os.remove(temp_html_file_path)
                        else: st.error("Failed to save/find temporary Pyvis HTML.")
                    except Exception as e_save: st.error(f"Error saving Pyvis graph: {e_save}"); st.session_state.pyvis_html_gn_gt = None
            elif st.session_state.static_string_image_url_gn_gt and not st.session_state.error_gn_gt: st.info("No detailed interaction data for Pyvis, static image might be available.")
            elif not st.session_state.error_gn_gt: st.warning("Could not retrieve network data/image. Check parameters.")
        else: st.warning("Please enter gene symbol and select organism.")
    if st.session_state.pyvis_html_gn_gt:
        st.markdown("---"); st.components.v1.html(st.session_state.pyvis_html_gn_gt, height=800, scrolling=True)
        st.download_button("⬇️ Download Interactive Network (HTML)", st.session_state.pyvis_html_gn_gt, f"{gene_symbol_gn}_interactive_network.html", "text/html", key="dl_pyvis_html_gt")
    if st.session_state.static_string_image_url_gn_gt:
        st.markdown("---"); st.subheader("Static Network Image (from STRING DB)")
        try: st.image(st.session_state.static_string_image_url_gn_gt, caption=f"Static STRING network for {gene_symbol_gn}", use_column_width=True)
        except Exception as e_img_display: st.warning(f"Could not display static image from URL: {e_img_display}")
        if st.button("Prepare Static PNG Download", key="prep_png_static_gn_gt"):
            with st.spinner("Downloading image data..."):
                try:
                    img_resp = requests.get(st.session_state.static_string_image_url_gn_gt, stream=True, timeout=15); img_resp.raise_for_status()
                    st.session_state.string_image_bytes_gn_gt = img_resp.content; st.session_state.string_image_filename_gn_gt = f"{gene_symbol_gn}_static_network.png"; st.success("Static PNG ready.")
                except requests.exceptions.RequestException as e_req: st.error(f"Failed to download static image: {e_req}")
                except Exception as e_gen: st.error(f"Unexpected error fetching static image: {e_gen}")
        if st.session_state.string_image_bytes_gn_gt: st.download_button("⬇️ Download Static Network (PNG)", st.session_state.string_image_bytes_gn_gt, st.session_state.string_image_filename_gn_gt, "image/png", key="dl_static_png_gt")
        if selected_organism_gn in ORGANISM_TO_TAXID and gene_symbol_gn:
            string_db_link = f"https://string-db.org/network/{ORGANISM_TO_TAXID[selected_organism_gn]}/{gene_symbol_gn.upper()}"; st.markdown(f"<a href='{string_db_link}' target='_blank' class='custom-link-button' style='display:inline-block; text-decoration:none; margin-top:10px;'>View on STRING DB</a>", unsafe_allow_html=True)

def render_blast_section():
    st.subheader("BLAST Search (NCBI)")
    st.markdown("Run a BLAST search against NCBI databases.")
    sequence = st.text_area("Enter FASTA sequence:", height=150, key="blast_sequence_gt")
    program = st.selectbox("BLAST Program:", ["blastn", "blastp", "blastx", "tblastn", "tblastx"], key="blast_program_gt")
    database = st.selectbox("Database:", ["nr", "nt", "refseq_protein", "pdb", "swissprot"], key="blast_database_gt")
    with st.expander("Advanced BLAST Options"):
        evalue_threshold = st.number_input("E-value threshold:", value=10.0, min_value=0.000001, max_value=1000.0, step=0.1, format="%.2e", key="blast_evalue_gt")
        max_hits = st.number_input("Max target sequences (hits):", value=50, min_value=1, max_value=500, step=10, key="blast_maxhits_gt")
    if st.button("Run BLAST", use_container_width=True, type="primary", key="blast_run_button_gt"):
        if not sequence.strip(): st.warning("Please enter a sequence."); return
        if not sequence.strip().startswith(">"): st.warning("Input sequence should be FASTA format. Adding default header."); sequence = ">QuerySequence\n" + sequence.strip()
        with st.spinner("Running BLAST (this may take a few minutes)..."):
            try:
                from Bio.Blast import NCBIWWW, NCBIXML; Entrez.email = "kaletejal05@mail.com"
                if Entrez.email == "kaletejal05@mail.com": st.sidebar.caption("INFO: Default Entrez email used for BLAST.")
                result_handle = NCBIWWW.qblast(program, database, sequence, expect=evalue_threshold, hitlist_size=max_hits)
                xml_output = result_handle.read(); result_handle.close()
                if not xml_output.strip(): st.error("BLAST query returned empty result."); return
                blast_records_iter = NCBIXML.parse(io.StringIO(xml_output)); results = []
                for blast_record in blast_records_iter:
                    if not blast_record.alignments: continue
                    for alignment in blast_record.alignments:
                        for hsp in alignment.hsps: results.append({"Hit ID": alignment.accession, "Description": alignment.hit_def, "Length": alignment.length, "Score": hsp.score, "Bits": hsp.bits, "E-value": hsp.expect, "Query Cover": f"{round((hsp.query_end - hsp.query_start + 1) / blast_record.query_length * 100, 1)}%", "Per. Ident": f"{round(hsp.identities / hsp.align_length * 100, 1)}%", "Gaps": hsp.gaps, "Alignment Length": hsp.align_length})
                if results:
                    df = pd.DataFrame(results); st.dataframe(df)
                    st.download_button("Download BLAST Results (CSV)", df.to_csv(index=False).encode('utf-8'), "blast_results.csv", "text/csv", key="dl_blast_csv_gt")
                    st.download_button("Download BLAST Results (XML)", xml_output.encode('utf-8'), "blast_results.xml", "application/xml", key="dl_blast_xml_gt")
                else: st.info("No significant hits found with current parameters.")
            except Exception as e: st.error(f"BLAST Error: {e}"); st.error("Check network, NCBI status, or query/parameters.")

def render_pathways_section():
    st.subheader("Pathway Analysis")
    st.markdown("<p style='font-size: 0.9em; color: #b0b0b0;'>Search for KEGG and Reactome pathways associated with a gene. Preview and download pathway images directly.</p>", unsafe_allow_html=True)
    path_col1, path_col2 = st.columns(2)
    with path_col1: gene_symbol_path = st.text_input("Enter Gene Symbol (e.g., TP53):", key="pathway_gene_symbol_gt").upper()
    with path_col2: organism_path = st.selectbox("Select Organism:", ORGANISM_OPTIONS, key="pathway_organism_gt", index=ORGANISM_OPTIONS.index("Homo sapiens") if "Homo sapiens" in ORGANISM_OPTIONS else 0)
    if 'kegg_pathways_data_gt' not in st.session_state: st.session_state.kegg_pathways_data_gt = None
    if 'kegg_error_data_gt' not in st.session_state: st.session_state.kegg_error_data_gt = None
    if 'reactome_pathways_data_gt' not in st.session_state: st.session_state.reactome_pathways_data_gt = None
    if 'reactome_error_data_gt' not in st.session_state: st.session_state.reactome_error_data_gt = None
    if 'current_pathway_gene_gt' not in st.session_state: st.session_state.current_pathway_gene_gt = ""
    if st.button("Fetch Pathways", key="fetch_pathways_btn_gt", use_container_width=True, type="primary"):
        if not gene_symbol_path: st.warning("Please enter a gene symbol.")
        else:
            st.session_state.current_pathway_gene_gt = gene_symbol_path
            with st.spinner(f"Fetching pathways for {gene_symbol_path} in {organism_path}..."):
                st.session_state.kegg_pathways_data_gt, st.session_state.kegg_error_data_gt = None, None
                st.session_state.reactome_pathways_data_gt, st.session_state.reactome_error_data_gt = None, None
                for key in list(st.session_state.keys()):
                    if key.startswith(("kegg_img_bytes_gt_", "kegg_img_err_gt_", "reactome_svg_data_gt_", "reactome_svg_err_gt_")): del st.session_state[key]
                kegg_data, kegg_err = fetch_kegg_pathways(gene_symbol_path, organism_path)
                st.session_state.kegg_pathways_data_gt = kegg_data; st.session_state.kegg_error_data_gt = kegg_err
                reactome_data, reactome_err = fetch_reactome_pathways(gene_symbol_path, organism_path)
                st.session_state.reactome_pathways_data_gt = reactome_data; st.session_state.reactome_error_data_gt = reactome_err

    if st.session_state.kegg_pathways_data_gt is not None or st.session_state.reactome_pathways_data_gt is not None or st.session_state.kegg_error_data_gt or st.session_state.reactome_error_data_gt:
        if st.session_state.get('kegg_error_data_gt'):
            st.error(f"KEGG Error: {st.session_state.kegg_error_data_gt}")
        if st.session_state.get('reactome_error_data_gt'):
            st.error(f"Reactome Error: {st.session_state.reactome_error_data_gt}")

    kegg_p_data = st.session_state.get('kegg_pathways_data_gt')
    if kegg_p_data and isinstance(kegg_p_data, list):
        st.subheader("KEGG Pathway Images")
        for p_idx, p_kegg in enumerate(kegg_p_data):
            if not all(k in p_kegg for k in ['Pathway ID', 'Name', 'URL']): continue
            st.markdown(f"**{p_kegg['Name']}** ([{p_kegg['Pathway ID']}]({p_kegg['URL']}))")
            img_bytes_key = f"kegg_img_bytes_gt_{p_kegg['Pathway ID']}"; img_err_key = f"kegg_img_err_gt_{p_kegg['Pathway ID']}"
            if st.button(f"Load Image for {p_kegg['Pathway ID']}", key=f"load_kegg_img_gt_{p_kegg['Pathway ID']}_{p_idx}", help="Click to load pathway diagram"):
                with st.spinner(f"Loading KEGG image for {p_kegg['Pathway ID']}..."):
                    img_bytes_kegg, img_err_kegg = fetch_kegg_pathway_image(p_kegg['Pathway ID'])
                    st.session_state[img_bytes_key] = img_bytes_kegg; st.session_state[img_err_key] = img_err_kegg
            if img_bytes_key in st.session_state:
                img_bytes_kegg_val = st.session_state[img_bytes_key]; img_err_kegg_val = st.session_state.get(img_err_key)
                if img_bytes_kegg_val:
                    try:
                        st.image(img_bytes_kegg_val, caption=f"KEGG Pathway: {p_kegg['Name']}", use_column_width=True)
                        st.download_button(f"⬇️ Download {p_kegg['Pathway ID']} Image", img_bytes_kegg_val, file_name=f"{p_kegg['Pathway ID']}_{st.session_state.current_pathway_gene_gt}.png", mime="image/png", key=f"dl_kegg_img_gt_{p_kegg['Pathway ID']}_{p_idx}_final")
                    except Exception as e_img_disp: st.error(f"Error displaying KEGG image {p_kegg['Pathway ID']}: {e_img_disp}")
                elif img_err_kegg_val: st.caption(f"Image for {p_kegg['Pathway ID']} not available: {img_err_kegg_val}")
            st.markdown("---")
    reactome_p_data = st.session_state.get('reactome_pathways_data_gt')
    if reactome_p_data and isinstance(reactome_p_data, list):
        st.subheader("Reactome Pathway Images (SVG)")
        for p_idx_r, p_reactome in enumerate(reactome_p_data):
            if not all(k in p_reactome for k in ['Pathway ID', 'Name', 'URL']): continue
            st.markdown(f"**{p_reactome['Name']}** ([{p_reactome['Pathway ID']}]({p_reactome['URL']}))")
            svg_url = f"https://reactome.org/ContentService/exporter/diagram/{p_reactome['Pathway ID']}.svg"
            st.markdown(f"[View SVG Diagram in New Tab]({svg_url})", unsafe_allow_html=True)
            svg_data_key = f"reactome_svg_data_gt_{p_reactome['Pathway ID']}"; svg_err_key = f"reactome_svg_err_gt_{p_reactome['Pathway ID']}"
            if st.button(f"Load SVG for {p_reactome['Pathway ID']}", key=f"load_reactome_svg_gt_{p_reactome['Pathway ID']}_{p_idx_r}", help="Click to load SVG diagram (can be slow)"):
                with st.spinner(f"Loading Reactome SVG for {p_reactome['Pathway ID']}..."):
                    try:
                        svg_response = requests.get(svg_url, timeout=20)
                        if svg_response.status_code == 200: st.session_state[svg_data_key] = svg_response.content.decode("utf-8"); st.session_state[svg_err_key] = None
                        else: st.session_state[svg_err_key] = f"SVG preview not available (Status: {svg_response.status_code})."; st.session_state[svg_data_key] = None
                    except requests.exceptions.Timeout: st.session_state[svg_err_key] = "Timeout when fetching SVG diagram."; st.session_state[svg_data_key] = None
                    except Exception as e_svg: st.session_state[svg_err_key] = f"SVG preview error: {e_svg}"; st.session_state[svg_data_key] = None
            if svg_data_key in st.session_state:
                svg_data_r_val = st.session_state[svg_data_key]; svg_err_r_val = st.session_state.get(svg_err_key)
                if svg_data_r_val:
                    st.markdown(f'<div style="border:1px solid #4CAF50; overflow:auto; background-color:white; max-height:500px; width:100%; padding:10px;">Preview for {p_reactome["Pathway ID"]}</div>', unsafe_allow_html=True)
                    st.download_button(f"⬇️ Download {p_reactome['Pathway ID']} SVG", svg_data_r_val, file_name=f"{p_reactome['Pathway ID']}_{st.session_state.current_pathway_gene_gt}.svg", mime="image/svg+xml", key=f"dl_reactome_svg_gt_{p_reactome['Pathway ID']}_{p_idx_r}_final")
                elif svg_err_r_val: st.caption(svg_err_r_val)
            st.markdown("---")

def render_primer_design_section():
    st.subheader("Primer Design (using Primer3)")
    st.markdown("Design PCR primers for a given DNA sequence.")
    if primer3 is None:
        st.error("Primer Design feature unavailable: `primer3-py` library or `primer3_core` not installed. Please install it (`pip install primer3-py`) and ensure primer3_core is in your PATH or accessible.")
        return

    template_sequence = st.text_area("Enter DNA Sequence Template (FASTA format or raw sequence):", height=150, key="primer_template_seq_gt", placeholder=">my_gene_id\nATGCATGCATGC...\nATGCATGCATGC...")
    st.markdown("##### Primer Design Parameters:")
    col1_pd, col2_pd, col3_pd = st.columns(3)
    with col1_pd:
        prod_min = st.number_input("Min Product Size (bp):", 100, 1000, 150, 10, key="primer_prod_min_gt")
        opt_size = st.number_input("Optimal Primer Size (nt):", 18, 25, 20, 1, key="primer_opt_size_gt")
        opt_tm = st.number_input("Optimal Primer Tm (°C):", 50.0, 70.0, 60.0, 0.5, format="%.1f", key="primer_opt_tm_gt")
    with col2_pd:
        prod_max = st.number_input("Max Product Size (bp):", prod_min + 50, 2000, 500, 10, key="primer_prod_max_gt")
        min_size = st.number_input("Min Primer Size (nt):", 15, opt_size, 18, 1, key="primer_min_size_gt")
        min_tm = st.number_input("Min Primer Tm (°C):", 45.0, opt_tm, 57.0, 0.5, format="%.1f", key="primer_min_tm_gt")
    with col3_pd:
        st.markdown("<div style='height:2.2em'></div>", unsafe_allow_html=True)
        max_size = st.number_input("Max Primer Size (nt):", opt_size, 30, 22, 1, key="primer_max_size_gt")
        max_tm = st.number_input("Max Primer Tm (°C):", opt_tm, 75.0, 63.0, 0.5, format="%.1f", key="primer_max_tm_gt")

    if st.button("Design Primers", key="design_primers_btn_gt", use_container_width=True, type="primary"):
        seq_id = "my_sequence"; actual_sequence = ""
        if not template_sequence.strip():
            st.warning("Please enter a DNA sequence.")
            return
        if template_sequence.startswith(">"):
            try:
                fasta_io = io.StringIO(template_sequence); record = next(SeqIO.parse(fasta_io, "fasta"))
                seq_id = record.id if record.id else "fasta_record"; actual_sequence = str(record.seq).upper()
                if not actual_sequence: st.error("FASTA record found but sequence is empty."); return
            except StopIteration: st.error("FASTA header found, but no sequence data followed."); return
            except Exception as e: st.error(f"Error parsing FASTA sequence: {e}."); return
        else:
            actual_sequence = template_sequence.strip().upper()

        if not actual_sequence:
            st.warning("Processed sequence is empty.")
            return

        param_errors = []
        if prod_min >= prod_max: param_errors.append("Min Product Size must be less than Max Product Size.")
        if min_size >= max_size: param_errors.append("Min Primer Size must be less than Max Primer Size.")
        if not (min_size <= opt_size <= max_size): param_errors.append("Optimal Primer Size must be within Min and Max Primer Size.")
        if min_tm >= max_tm: param_errors.append("Min Primer Tm must be less than Max Primer Tm.")
        if not (min_tm <= opt_tm <= max_tm): param_errors.append("Optimal Primer Tm must be within Min and Max Primer Tm.")
        if param_errors:
            for err_msg in param_errors: st.error(err_msg)
            return

        with st.spinner("Designing primers with Primer3..."):
            results, error = design_primers_with_primer3(
                seq_id, actual_sequence,
                prod_min, prod_max,
                opt_size, min_size, max_size,
                opt_tm, min_tm, max_tm
            )

        if error and not results:
            st.error(f"Primer Design Error: {error}")
        elif results:
            num_pairs_returned = results.get('PRIMER_PAIR_NUM_RETURNED', 0)
            if error and num_pairs_returned == 0 :
                 st.info(f"Note on Primer Design: {error}")
            elif not error and num_pairs_returned == 0:
                 st.info("No primer pairs found with current parameters. The sequence might be unsuitable (e.g., too short, high N content, repetitive) or parameters too restrictive.")

            if num_pairs_returned > 0:
                st.success(f"Primer design successful! Found {num_pairs_returned} primer pair(s).")
                primer_data = []
                for i in range(num_pairs_returned):
                    pair_info = {
                        "Pair #": i + 1,
                        "Fwd Sequence": results.get(f'PRIMER_LEFT_{i}_SEQUENCE', 'N/A'),
                        "Fwd Start": results.get(f'PRIMER_LEFT_{i}', (None, None))[0],
                        "Fwd Length (nt)": results.get(f'PRIMER_LEFT_{i}', (None, None))[1],
                        "Fwd Tm (°C)": results.get(f'PRIMER_LEFT_{i}_TM', 'N/A'),
                        "Fwd GC (%)": results.get(f'PRIMER_LEFT_{i}_GC_PERCENT', 'N/A'),
                        "Rev Sequence": results.get(f'PRIMER_RIGHT_{i}_SEQUENCE', 'N/A'),
                        "Rev Start": results.get(f'PRIMER_RIGHT_{i}', (None, None))[0],
                        "Rev Length (nt)": results.get(f'PRIMER_RIGHT_{i}', (None, None))[1],
                        "Rev Tm (°C)": results.get(f'PRIMER_RIGHT_{i}_TM', 'N/A'),
                        "Rev GC (%)": results.get(f'PRIMER_RIGHT_{i}_GC_PERCENT', 'N/A'),
                        "Product Size (bp)": results.get(f'PRIMER_PAIR_{i}_PRODUCT_SIZE', 'N/A'),
                    }
                    primer_data.append(pair_info)
                df_primers = pd.DataFrame(primer_data)
                st.dataframe(df_primers.style.format({
                    "Fwd Tm (°C)": "{:.2f}", "Fwd GC (%)": "{:.2f}",
                    "Rev Tm (°C)": "{:.2f}", "Rev GC (%)": "{:.2f}"
                }))
                csv = df_primers.to_csv(index=False).encode('utf-8')
                st.download_button("Download Primer Results (CSV)", csv, f"primers_{seq_id}.csv", "text/csv", key=f'download-primers-csv-{seq_id}')

            with st.expander("View Full Primer3 Output (Raw Dictionary)"):
                st.json(results)
        else:
             st.error("Primer design failed unexpectedly. No results or error message returned from the design function.")


# --- New Page for Genomic Toolkit ---
def render_genomic_toolkit_page():
    st.markdown("<h1 class='main-header'>Genomic Toolkit</h1>", unsafe_allow_html=True)
    st.markdown("<p style='font-size: 1.1em; color: #e0e0e0;'>A suite of tools for comprehensive gene and genome analysis.</p>", unsafe_allow_html=True)
    st.markdown("---")
    tool_tabs = st.tabs([
        " Gene Search", " Gene Network", " Pathways",
        " BLAST Search", " Primer Design"
    ])
    with tool_tabs[0]: render_gene_search_section()
    with tool_tabs[1]: render_gene_network_section()
    with tool_tabs[2]: render_pathways_section()
    with tool_tabs[3]: render_blast_section()
    with tool_tabs[4]: render_primer_design_section()

# --- FASTA Analysis Tools Page ---
def render_fasta_analysis_tools_page():
    st.markdown("<h1 class='main-header'>FASTA Analysis Utilities</h1>", unsafe_allow_html=True)
    st.markdown("<p style='font-size: 1.1em; color: #e0e0e0;'>Tools for analyzing FASTA sequences and related data.</p>", unsafe_allow_html=True)
    st.markdown("---")

    st.subheader(" GC Profile Visualizer")
    fasta_input_gc = st.text_area("Paste FASTA sequence(s) for GC Profile:", height=150, key="gc_profile_input_fasta", help="One or more DNA/RNA sequences in FASTA format.")
    uploaded_file_gc = st.file_uploader("Or upload a FASTA file for GC Profile:", type=['fasta', 'fa', 'fna', 'rna'], key="gc_profile_upload_fasta")
    window_size_gc = st.number_input("Window Size (bp/nt):", min_value=10, value=50, step=10, key="gc_window_size_input", help="e.g., 50-200 bp/nt.")
    if st.button("Generate GC Profile Plot", key="gc_profile_btn_generate", use_container_width=True, type="primary"):
        sequence_data_gc_str = None
        if uploaded_file_gc:
            try: sequence_data_gc_str = uploaded_file_gc.read().decode()
            except Exception as e: st.error(f"Error reading file: {e}")
        elif fasta_input_gc.strip():
            sequence_data_gc_str = fasta_input_gc
        if sequence_data_gc_str:
            parsed_sequences = parse_fasta_sequences_for_gc(sequence_data_gc_str)
            if not parsed_sequences:
                st.error("No valid DNA/RNA sequences found. Check FASTA format/content.")
            else:
                num_to_plot = len(parsed_sequences)
                if num_to_plot > 5:
                    st.warning(f"Found {num_to_plot} sequences. Plotting first 5 for performance.")
                    parsed_sequences = parsed_sequences[:5]
                for i, seq_info in enumerate(parsed_sequences):
                    st.markdown(f"#### GC Profile: `{seq_info['id']}` (Length: {len(seq_info['seq'])} bp/nt)")
                    if len(seq_info['seq']) < window_size_gc:
                        st.warning(f"Sequence '{seq_info['id']}' ({len(seq_info['seq'])} bp/nt) is shorter than window size ({window_size_gc} bp/nt). Cannot plot."); continue
                    with st.spinner(f"Calculating for {seq_info['id']}..."):
                        positions, gc_values = calculate_gc_profile(seq_info['seq'], window_size_gc)
                    if not positions:
                        st.warning(f"Could not calculate GC profile for '{seq_info['id']}'. Ensure it contains standard DNA/RNA bases."); continue
                    fig = go.Figure(data=[go.Scatter(x=positions, y=gc_values, mode='lines', name='GC %', line=dict(color='#4CAF50', width=2))])
                    fig.update_layout(title=f"GC Content: {seq_info['id']}", xaxis_title="Position (bp/nt midpoint)", yaxis_title="GC Content (%)", yaxis_range=[0,100], plot_bgcolor="#23272f", paper_bgcolor="#23272f", font=dict(color="#eaffea"), height=400)
                    st.plotly_chart(fig, use_container_width=True)
                    if i < len(parsed_sequences) - 1: st.markdown("<hr style='border-top: 1px dashed #00bcd4;'>", unsafe_allow_html=True)
        else: st.warning("Please paste a sequence or upload a FASTA file.")
    st.markdown("---")

    st.subheader(" Base Composition Analyzer")
    dna_sequence_comp_input = st.text_area("Paste single DNA/RNA sequence for Composition:", height=100, key="base_comp_input_dna", help="Single raw DNA/RNA sequence (not FASTA).")
    if st.button("Analyze Base Composition", key="base_comp_btn_analyze", use_container_width=True, type="primary"):
        if dna_sequence_comp_input.strip():
            sequence_to_analyze = dna_sequence_comp_input.strip()
            counts, total_length = calculate_nucleotide_composition(sequence_to_analyze)
            if total_length == 0: st.warning("Input sequence is empty.")
            else:
                labels = [base for base, count in counts.items() if count > 0]
                values = [count for base, count in counts.items() if count > 0]
                if not values: st.warning("No standard bases (A, T, G, C, U) or 'Other' found to plot.")
                else:
                    fig_pie = px.pie(names=labels, values=values, title=f"Nucleotide Composition (Total length: {total_length} bp/nt)",
                                     color_discrete_map={'A':'#63A074', 'T':'#D6A163', 'G':'#6B8CD6', 'C':'#B76E6E', 'U':'#8E63D6', 'Other':'#9E9E9E'})
                    fig_pie.update_layout(plot_bgcolor="#23272f", paper_bgcolor="#23272f", font=dict(color="#eaffea"), legend_title_text='Bases')
                    st.plotly_chart(fig_pie, use_container_width=True)
                    st.markdown("##### Base Counts:")
                    count_str = []
                    for base, count in counts.items():
                        if count > 0:
                            percentage = (count / total_length) * 100 if total_length > 0 else 0
                            count_str.append(f"**{base}**: {count} ({percentage:.2f}%)")
                    st.markdown("      ".join(count_str), unsafe_allow_html=True)
        else: st.warning("Please paste a DNA/RNA sequence.")
    st.markdown("---")

    st.subheader(" ORF Finder")
    st.markdown("Find Open Reading Frames (ORFs) in DNA sequences.")
    orf_input_type = st.radio(
        "Input type for ORF Finder:",
        ("Upload FASTA file", "Paste raw DNA sequence"),
        key="orf_input_type_fasta_utils"
    )
    orf_sequence_input_str = ""
    uploaded_orf_file = None
    seq_id_for_orf = "input_sequence"
    if orf_input_type == "Upload FASTA file":
        uploaded_orf_file = st.file_uploader(
            "Upload FASTA file for ORF Finding (uses first sequence):",
            type=['fasta', 'fa', 'fna'],
            key="orf_upload_fasta"
        )
        if uploaded_orf_file:
            try:
                fasta_record = next(SeqIO.parse(io.StringIO(uploaded_orf_file.read().decode()), "fasta"))
                orf_sequence_input_str = str(fasta_record.seq)
                seq_id_for_orf = fasta_record.id if fasta_record.id else "fasta_record"
                st.info(f"Using sequence ID: {seq_id_for_orf} (Length: {len(orf_sequence_input_str)} bp)")
            except StopIteration:
                st.error("No sequences found in the uploaded FASTA file.")
                orf_sequence_input_str = ""
            except Exception as e:
                st.error(f"Error reading FASTA file: {e}")
                orf_sequence_input_str = ""
    else:
        orf_sequence_input_str_raw = st.text_area(
            "Paste raw DNA sequence (no headers):",
            height=150,
            key="orf_raw_dna_input"
        )
        if orf_sequence_input_str_raw:
            orf_sequence_input_str = orf_sequence_input_str_raw.strip().upper()
            if orf_sequence_input_str:
                 st.info(f"Processing pasted sequence (Length: {len(orf_sequence_input_str)} bp)")

    min_orf_length_aa = st.number_input(
        "Minimum ORF length (amino acids):",
        min_value=10, value=30, step=5, key="orf_min_len_aa_v2"
    )
    min_orf_length_bp = (min_orf_length_aa * 3) + 3

    if st.button("Find ORFs", key="orf_btn_find_v2", use_container_width=True, type="primary"):
        if not orf_sequence_input_str:
            st.warning("Please provide a DNA sequence.")
        else:
            valid_dna_chars = "ATGC"
            cleaned_sequence = "".join([char for char in orf_sequence_input_str.upper() if char in valid_dna_chars])
            if len(cleaned_sequence) != len(orf_sequence_input_str):
                st.warning(f"Non-DNA characters were removed. Original length: {len(orf_sequence_input_str)}, Processed length: {len(cleaned_sequence)}")
            if not cleaned_sequence:
                st.error("Sequence contains no valid DNA bases (A, T, G, C) after cleaning.")
            elif len(cleaned_sequence) < min_orf_length_bp:
                 st.warning(f"Cleaned sequence length ({len(cleaned_sequence)} bp) is less than the minimum ORF length required ({min_orf_length_bp} bp for {min_orf_length_aa} aa).")
            else:
                with st.spinner("Finding ORFs..."):
                    found_orfs_list = find_orfs(cleaned_sequence, min_orf_len_bp=min_orf_length_bp)
                if found_orfs_list:
                    st.success(f"Found {len(found_orfs_list)} ORFs meeting the criteria.")
                    df_orfs = pd.DataFrame(found_orfs_list)
                    df_orfs_display = df_orfs.copy()
                    df_orfs_display["Protein Preview"] = df_orfs_display["Protein Sequence"].apply(lambda x: x[:30] + '...' if len(x) > 30 else x)
                    st.dataframe(df_orfs_display[[
                        "Frame", "Start (0-indexed)", "End (0-indexed)",
                        "Length (bp)", "Protein Length (aa)", "Protein Preview"
                    ]], height=300)
                    dna_fasta_orfs = ""
                    for i, orf in enumerate(found_orfs_list):
                        dna_fasta_orfs += f">{seq_id_for_orf}_ORF_{i+1}_Frame{orf['Frame']}_Start{orf['Start (0-indexed)']}_Len{orf['Length (bp)']}\n{orf['DNA Sequence']}\n"
                    st.download_button(
                        label="⬇️ Download ORF DNA Sequences (FASTA)",
                        data=dna_fasta_orfs,
                        file_name=f"{seq_id_for_orf}_orfs_dna.fasta",
                        mime="text/fasta",
                        key="download_orf_dna_fasta_v2"
                    )
                    protein_fasta_orfs = ""
                    for i, orf in enumerate(found_orfs_list):
                        if not orf['Protein Sequence'].startswith("[Translation Error"):
                            protein_fasta_orfs += f">{seq_id_for_orf}_ORF_{i+1}_Frame{orf['Frame']}_Start{orf['Start (0-indexed)']}_ProtLen{orf['Protein Length (aa)']}\n{orf['Protein Sequence']}\n"
                    if protein_fasta_orfs:
                        st.download_button(
                            label="⬇️ Download ORF Protein Sequences (FASTA)",
                            data=protein_fasta_orfs,
                            file_name=f"{seq_id_for_orf}_orfs_protein.fasta",
                            mime="text/fasta",
                            key="download_orf_protein_fasta_v2"
                        )
                else:
                    st.info("No ORFs found matching the specified criteria for the provided sequence.")
    st.markdown("---")

    st.subheader(" 3D Structure Viewer (PDB)")
    st.markdown("Visualize 3D protein structures from PDB IDs.")
    if py3Dmol is None:
        st.error("`py3Dmol` library is not available. This feature is disabled. Please install it: `pip install py3Dmol`")
    else:
        pdb_id_input = st.text_input("Enter PDB ID (e.g., 1MOL, 6LU7):", key="pdb_id_3d_viewer").strip().upper()
        viewer_height = st.slider("Viewer Height (pixels):", 300, 800, 500, 50, key="3d_viewer_height")
        style_options = ["cartoon", "stick", "sphere", "line", "cross"]
        selected_style = st.selectbox("Select Structure Style:", style_options, key="3d_style_select")
        color_options = {
            "Chain": "chain", "Spectrum": "spectrum", "SS (Secondary Structure)": "ssJmol",
            "Residue Type": "residue", "Default (White)": "white"
        }
        selected_color_key = st.selectbox("Select Color Scheme:", list(color_options.keys()), key="3d_color_select")
        actual_color_scheme = color_options[selected_color_key]

        if st.button("Visualize Structure", key="btn_visualize_3d", use_container_width=True, type="primary"):
            if not pdb_id_input:
                st.warning("Please enter a PDB ID.")
            elif len(pdb_id_input) != 4:
                st.warning("PDB ID should typically be 4 characters long.")
            else:
                with st.spinner(f"Fetching and rendering PDB ID: {pdb_id_input}..."):
                    try:
                        view = py3Dmol.view(width="100%", height=viewer_height, viewergrid=(1,1), js="https://3Dmol.org/build/3Dmol-min.js")
                        view.addModel(requests.get(f"https://files.rcsb.org/download/{pdb_id_input}.pdb").text, "pdb")
                        
                        style_spec = {selected_style: {}}
                        if actual_color_scheme == 'ssJmol' and selected_style == 'cartoon':
                            style_spec[selected_style]['colorscheme'] = 'ssJmol'
                        elif actual_color_scheme == 'white':
                            style_spec[selected_style]['color'] = 'white'
                        else: 
                            style_spec[selected_style]['colorscheme'] = f'{actual_color_scheme}Carbon' if selected_style != 'cartoon' else actual_color_scheme
                            if selected_style == 'cartoon' and actual_color_scheme not in ['chain', 'spectrum', 'residue', 'ssJmol']: 
                                style_spec[selected_style]['color'] = actual_color_scheme 
                            elif selected_style == 'cartoon' and actual_color_scheme in ['chain', 'spectrum', 'residue']:
                                style_spec[selected_style]['colorscheme'] = actual_color_scheme


                        view.setStyle({}, style_spec)
                        view.setBackgroundColor('0x23272f') 
                        view.zoomTo()
                        view.show() 
                        st.components.v1.html(view._make_html(), height=viewer_height + 20, scrolling=False)
                        st.caption(f"Displaying PDB: {pdb_id_input}. Rotate with mouse, scroll to zoom.")
                    except requests.exceptions.RequestException as e:
                        st.error(f"Failed to fetch PDB data for {pdb_id_input}. Error: {e}. Check PDB ID and internet connection.")
                    except AttributeError as e:
                        st.error(f"Error with py3Dmol rendering: {e}. The library might have updated or there's an issue with its methods.")
                    except Exception as e:
                        st.error(f"An unexpected error occurred: {e}")
                        st.exception(e)
    st.markdown("---")

    st.subheader(" Protein Residue Counter")
    protein_sequence_aa = st.text_area("Paste single protein sequence for AA Frequency:", height=100, key="aa_freq_input_protein", help="Single raw protein sequence (not FASTA).")
    if st.button("Calculate Amino Acid Frequencies", key="aa_freq_btn_calculate", use_container_width=True, type="primary"):
        protein_sequence_aa_stripped = protein_sequence_aa.strip()
        if protein_sequence_aa_stripped:
            with st.spinner("Calculating amino acid frequencies..."):
                aa_composition_df = calculate_amino_acid_composition(protein_sequence_aa_stripped)
            if not aa_composition_df.empty:
                st.markdown("##### Amino Acid Composition:")
                plot_df = aa_composition_df[aa_composition_df['Percentage'] > 0.001];
                if plot_df.empty and not aa_composition_df[aa_composition_df['Count'] > 0].empty: plot_df = aa_composition_df[aa_composition_df['Count'] > 0]
                if plot_df.empty and len(protein_sequence_aa_stripped) > 0 : st.info(f"Sequence of length {len(protein_sequence_aa_stripped)} had no standard amino acids.")
                elif not plot_df.empty:
                    fig_aa_bar = px.bar(plot_df.sort_values(by="Percentage", ascending=False), x="Amino Acid", y="Percentage", title=f"Amino Acid Frequency (Total residues: {len(protein_sequence_aa_stripped)})", labels={"Percentage": "Frequency (%)", "Amino Acid": "Amino Acid"}, color="Amino Acid", text_auto='.2f')
                    fig_aa_bar.update_layout(plot_bgcolor="#23272f", paper_bgcolor="#23272f", font=dict(color="#eaffea"), yaxis_title="Frequency (%)", showlegend=False)
                    fig_aa_bar.update_traces(texttemplate='%{y:.2f}%', textposition='outside'); st.plotly_chart(fig_aa_bar, use_container_width=True)
                st.dataframe(aa_composition_df.style.format({"Percentage": "{:.2f}%"}), column_config={"Amino Acid": st.column_config.TextColumn("Amino Acid", width="small"), "Count": st.column_config.NumberColumn("Count", width="small"), "Percentage": st.column_config.NumberColumn("Percentage", format="%.2f%%", width="medium")}, use_container_width=True, hide_index=True)
        else: st.warning("Please paste a protein sequence.")
    st.markdown("---")

    st.subheader(" Expression Pattern Mapper")
    uploaded_csv_heatmap = st.file_uploader("Upload CSV for Expression Heatmap:", type=['csv'], key="heatmap_csv_upload_expr")
    if uploaded_csv_heatmap is not None:
        try:
            df_original = pd.read_csv(uploaded_csv_heatmap)
            if df_original.empty: st.warning("Uploaded CSV is empty.")
            else:
                st.markdown("##### Heatmap Options")
                all_columns = df_original.columns.tolist(); default_index_col_guess = None
                if all_columns:
                    if pd.api.types.is_string_dtype(df_original[all_columns[0]]): default_index_col_guess = all_columns[0]
                    else:
                        for col_name_idx in all_columns:
                            if pd.api.types.is_string_dtype(df_original[col_name_idx]): default_index_col_guess = col_name_idx; break
                        if default_index_col_guess is None: default_index_col_guess = all_columns[0]
                index_col = st.selectbox("Select Index Column (Genes/Features):", all_columns, index=all_columns.index(default_index_col_guess) if default_index_col_guess and default_index_col_guess in all_columns else 0, key="heatmap_index_col")
                data_cols_options = [col for col in all_columns if col != index_col]
                if not data_cols_options and index_col: st.warning("No data columns available after selecting index. CSV needs >=2 columns.")
                selected_data_cols = st.multiselect("Select Data Columns (Samples/Conditions):", data_cols_options, default=data_cols_options, key="heatmap_data_cols")
                available_colorscales = px.colors.named_colorscales(); default_colorscale_name = "viridis" if "viridis" in available_colorscales else available_colorscales[0]
                try: default_colorscale_idx = available_colorscales.index(default_colorscale_name)
                except ValueError: default_colorscale_idx = 0
                colorscale = st.selectbox("Select Colorscale:", available_colorscales, index=default_colorscale_idx, key="heatmap_colorscale")
                col_norm1, col_norm2 = st.columns(2)
                with col_norm1: normalize_zscore_gene = st.checkbox("Normalize by Gene (Z-score rows)", key="heatmap_norm_gene_v2")
                with col_norm2: normalize_zscore_sample = st.checkbox("Normalize by Sample (Z-score columns)", key="heatmap_norm_sample_v2", disabled=st.session_state.get("heatmap_norm_gene_v2", False))
                if st.session_state.get("heatmap_norm_gene_v2", False) and st.session_state.get("heatmap_norm_sample_v2", False): st.session_state.heatmap_norm_sample_v2 = False
                if st.button("Generate Expression Heatmap", key="heatmap_btn_generate_expr_final", use_container_width=True, type="primary"):
                    if not selected_data_cols: st.error("Please select at least one data column for heatmap.")
                    else:
                        with st.spinner("Generating heatmap..."):
                            try:
                                df_processed = df_original.copy()
                                if index_col and index_col in df_processed.columns:
                                     if df_processed[index_col].duplicated().any(): st.warning(f"Duplicate values in index column ('{index_col}').")
                                     df_processed = df_processed.set_index(index_col)
                                else: st.error(f"Index column '{index_col}' not found."); return
                                df_heatmap_data = df_processed[selected_data_cols].copy()
                                failed_numeric_conversion = []
                                for col in df_heatmap_data.columns:
                                    original_dtype = df_heatmap_data[col].dtype; df_heatmap_data[col] = pd.to_numeric(df_heatmap_data[col], errors='coerce')
                                    if df_heatmap_data[col].isnull().all() and not pd.api.types.is_numeric_dtype(original_dtype): failed_numeric_conversion.append(col)
                                if failed_numeric_conversion: st.warning(f"Columns excluded (all NaN): {', '.join(failed_numeric_conversion)}"); df_heatmap_data = df_heatmap_data.drop(columns=failed_numeric_conversion, errors='ignore')
                                df_heatmap_data.dropna(how='all', axis=0, inplace=True); df_heatmap_data.dropna(how='all', axis=1, inplace=True)
                                if df_heatmap_data.empty: st.error("No valid numeric data to plot."); return
                                do_norm_gene = st.session_state.heatmap_norm_gene_v2; do_norm_sample = st.session_state.heatmap_norm_sample_v2
                                if do_norm_gene: df_heatmap_data = df_heatmap_data.apply(lambda x: zscore(x, nan_policy='omit') if x.notna().sum() > 1 else x, axis=1).fillna(0)
                                elif do_norm_sample: df_heatmap_data = df_heatmap_data.apply(lambda x: zscore(x, nan_policy='omit') if x.notna().sum() > 1 else x, axis=0).fillna(0)
                                fig_heatmap = go.Figure(data=go.Heatmap(z=df_heatmap_data.values, x=df_heatmap_data.columns.tolist(), y=df_heatmap_data.index.tolist(), colorscale=colorscale, colorbar_title="Value", hoverongaps=False ))
                                title_suffix = " (Z-score by Gene)" if do_norm_gene else " (Z-score by Sample)" if do_norm_sample else ""
                                fig_heatmap.update_layout(title=f"Expression Heatmap{title_suffix}", xaxis_title="Samples / Conditions", yaxis_title="Genes / Features", plot_bgcolor="#23272f", paper_bgcolor="#23272f", font=dict(color="#eaffea"), height=max(400, len(df_heatmap_data.index) * 18 + 150), xaxis=dict(tickangle=-45) )
                                st.plotly_chart(fig_heatmap, use_container_width=True)
                                with st.expander("View Processed Data Table for Heatmap"):
                                    table_cmap_map = {'reds': 'Reds', 'ylorrd': 'Reds', 'blues': 'Blues', 'ylgnbu': 'Blues', 'greens': 'Greens', 'coolwarm': 'coolwarm', 'plasma': 'plasma', 'magma': 'magma'}
                                    table_cmap = table_cmap_map.get(colorscale.lower(), 'viridis')
                                    st.dataframe(df_heatmap_data.style.background_gradient(cmap=table_cmap, axis=None).format("{:.3f}"))
                            except Exception as e_process: st.error(f"Error generating heatmap: {e_process}"); st.exception(e_process)
        except pd.errors.EmptyDataError: st.error("Uploaded CSV is empty or unparsable.")
        except Exception as e_read: st.error(f"Error reading/processing CSV: {e_read}"); st.exception(e_read)
    st.markdown("---")

    st.subheader(" Codon Usage Analyzer")
    st.markdown("Calculate codon usage from coding DNA sequences (CDS). Input can be a FASTA file or raw sequences.")
    codon_usage_input_type = st.radio("Input type for Codon Usage:", ("Upload FASTA file of CDS", "Paste raw CDS sequences (one per line)"), key="codon_usage_input_type_gt")
    cds_sequences_str_list = []
    if codon_usage_input_type == "Upload FASTA file of CDS":
        uploaded_cds_file = st.file_uploader("Upload FASTA file containing CDS:", type=['fasta', 'fa', 'fna'], key="codon_usage_upload_cds")
        if uploaded_cds_file:
            try:
                fasta_content = uploaded_cds_file.read().decode()
                fasta_io = io.StringIO(fasta_content)
                for record in SeqIO.parse(fasta_io, "fasta"):
                    if record.seq: cds_sequences_str_list.append(str(record.seq))
                if not cds_sequences_str_list and fasta_content: st.warning("FASTA parsed, but no sequences extracted.")
            except Exception as e: st.error(f"Error reading/parsing FASTA: {e}")
    else:
        raw_cds_input = st.text_area("Paste CDS sequences (one per line, no headers):", height=150, key="codon_usage_raw_cds_input", placeholder="ATGCGTAGCATGC...\nGGCTAGCTAGCT...")
        if raw_cds_input.strip(): cds_sequences_str_list = [seq.strip() for seq in raw_cds_input.splitlines() if seq.strip()]
    if st.button("Calculate Codon Usage", key="codon_usage_btn_calculate", use_container_width=True, type="primary"):
        if not cds_sequences_str_list: st.warning("Please provide CDS sequences.")
        else:
            with st.spinner("Calculating codon usage..."):
                df_usage, aa_codon_counts, error_msg = calculate_codon_usage_from_cds(cds_sequences_str_list)
            if error_msg: st.warning(f"Calculation Notes: {error_msg}")
            if not df_usage.empty:
                st.success(f"Codon usage calculated for {len(cds_sequences_str_list)} valid sequence(s) submitted.")
                st.markdown("##### Codon Usage Table:"); st.dataframe(df_usage.style.format({"Frequency (‰)": "{:.2f}"}), use_container_width=True)
                st.markdown("##### Amino Acid - Codon Preferences (Counts):")
                aa_to_plot = [aa for aa, codons in aa_codon_counts.items() if codons and aa != "?"]
                if aa_to_plot:
                    num_cols_aa_plot = st.slider("Columns for AA plots:", 2, 5, 3, key="aa_plot_cols")
                    plot_cols = st.columns(num_cols_aa_plot)
                    for i, aa in enumerate(sorted(aa_to_plot)):
                        with plot_cols[i % num_cols_aa_plot]:
                            codon_data = aa_codon_counts[aa]
                            if codon_data:
                                df_aa_codons = pd.DataFrame(list(codon_data.items()), columns=['Codon', 'Count']).sort_values(by="Count", ascending=False)
                                fig_aa = px.bar(df_aa_codons, x='Codon', y='Count', title=f"'{aa}' Codons")
                                fig_aa.update_layout(height=300, margin=dict(l=20, r=20, t=40, b=20), plot_bgcolor="#23272f", paper_bgcolor="#23272f", font=dict(color="#eaffea"))
                                st.plotly_chart(fig_aa, use_container_width=True)
                else: st.info("No specific amino acid codon counts to plot.")
            elif not error_msg: st.error("Could not calculate codon usage. Ensure valid CDS.")
    st.markdown("---")

    st.subheader(" Phylogenetic Tree Viewer (Newick)")
    st.markdown("Upload a phylogenetic tree in Newick format to visualize it.")
    uploaded_newick_file = st.file_uploader("Upload Newick tree file (.nwk, .newick, .tree, .txt):", type=['nwk', 'newick', 'tree', 'txt'], key="phylo_tree_upload_newick")
    if uploaded_newick_file:
        try:
            newick_string = uploaded_newick_file.read().decode().strip()
            if not newick_string: st.warning("Uploaded Newick file is empty.")
            else:
                with st.spinner("Parsing and rendering tree..."):
                    tree, error_parse = parse_newick_tree(newick_string)
                if error_parse: st.error(error_parse)
                elif tree:
                    st.success("Tree parsed successfully!")
                    try:
                        fig_plotly_tree = plot_phylo_tree_plotly(tree)
                        if fig_plotly_tree: st.plotly_chart(fig_plotly_tree, use_container_width=True)
                        else: st.warning("Could not generate Plotly tree visualization.")
                    except Exception as e_plotly:
                        st.error(f"Error generating Plotly tree: {e_plotly}"); st.exception(e_plotly)
                        st.info("Plotly visualization might fail for complex trees. Consider an external viewer.")
                    with st.expander("View Newick String"):
                        st.text_area("Newick Content:", newick_string, height=100, disabled=True)
                else: st.error("Failed to parse or render tree.")
        except Exception as e_file: st.error(f"Error reading uploaded file: {e_file}")
    st.markdown("---")


# --- Metagenomics Tools Page ---
def render_metagenomics_tools_page():
    st.markdown("<h1 class='main-header'>Metagenomics Resources</h1>", unsafe_allow_html=True)
    st.markdown("<p style='font-size: 1.1em; color: #e0e0e0;'>Explore curated databases and tools for metagenomics research. Click 'Open Tool' to attempt embedding or visit the external site.</p>", unsafe_allow_html=True)
    if 'metagen_tool_iframe_url' not in st.session_state: st.session_state.metagen_tool_iframe_url = None
    if 'metagen_tool_iframe_name' not in st.session_state: st.session_state.metagen_tool_iframe_name = None
    tool_categories = {
        "General Analysis Platforms & Databases": [{"name": "MG-RAST", "description": "A platform for annotation of metagenomic data...", "url": "https://www.mg-rast.org/", "iframe_able": True }, {"name": "MGnify", "description": "MGnify is a freely available hub for the analysis...", "url": "https://www.ebi.ac.uk/metagenomics/", "iframe_able": True}, {"name": "IMG/M", "description": "The Integrated Microbial Genomes (IMG) system...", "url": "https://img.jgi.doe.gov/m/", "iframe_able": True }, {"name": "JGI GOLD", "description": "Genomes Online Database, providing access to metadata.", "url": "https://gold.jgi.doe.gov/", "iframe_able": True}],
        "Taxonomic Profiling & Classification": [{"name": "MetaPhlAn", "description": "Computational tool for profiling microbial communities...", "url": "https://huttenhower.sph.harvard.edu/metaphlan/", "iframe_able": False }, {"name": "Kraken2 / Bracken", "description": "Kraken2 for fast taxonomic classification...", "url": "https://ccb.jhu.edu/software/kraken2/", "iframe_able": False}, {"name": "Kaiju", "description": "Fast and sensitive taxonomic classification...", "url": "http://kaiju.binf.ku.dk/", "iframe_able": True}],
        "Functional Annotation Tools": [{"name": "MetaGeneMark", "description": "Tool for identifying protein-coding regions...", "url": "http://topaz.gatech.edu/GeneMark/metagenome/", "iframe_able": True }, {"name": "eggNOG-mapper", "description": "Fast functional annotation of sequences...", "url": "http://eggnog-mapper.embl.de/", "iframe_able": True}, {"name": "KofamKOALA", "description": "KEGG Orthology (KO) assignment based on HMM profiles...", "url": "https://www.genome.jp/tools/kofamkoala/", "iframe_able": True}, {"name": "dbCAN", "description": "Automated carbohydrate-active enzyme (CAZyme) annotation...", "url": "http://bcb.unl.edu/dbCAN2/", "iframe_able": True}],
        "Metagenome Assembly & Binning Tools (Command-line)": [{"name": "MEGAHIT", "description": "A very fast and memory-efficient metagenome assembler...", "url": "https://github.com/voutcn/megahit", "iframe_able": False}, {"name": "metaSPAdes", "description": "A metagenomic assembler designed for complex datasets...", "url": "http://cab.spbu.ru/software/spades/#meta", "iframe_able": False}, {"name": "MetaBAT2", "description": "An efficient tool for reconstructing genomes...", "url": "https://bitbucket.org/berkeleylab/metabat2/src/master/", "iframe_able": False}],
        "Specialized Databases & Resources": [{"name": "GMGC (Global Microbial Gene Catalog)", "description": "A comprehensive catalog of genes from metagenomes...", "url": "https://gmgc.embl.de/", "iframe_able": True }, {"name": "Human Microbiome Project (HMP)", "description": "NIH's Human Microbiome Project provides microbiota data...", "url": "https://hmpdacc.org/", "iframe_able": True}, {"name": "iHMP (Integrated HMP)", "description": "Longitudinal multi-omics studies of microbiome dynamics...", "url": "https://www.hmpdacc.org/ihmp/", "iframe_able": True}, {"name": "CARD (Antibiotic Resistance)", "description": "Comprehensive Antibiotic Resistance Database...", "url": "https://card.mcmaster.ca/", "iframe_able": True}, {"name": "ResFinder", "description": "Web tool for identifying acquired antimicrobial resistance genes...", "url": "https://cge.food.dtu.dk/services/ResFinder/", "iframe_able": True}, {"name": "NMPFamsDB", "description": "Database of metagenome- and metatranscriptome-derived protein families...", "url": "http://nmpfams.org/", "iframe_able": True }]
    }
    def display_tool_section_cards(title, tools_list):
        st.markdown(f"<h3 class='section-subheader'>{title}</h3>", unsafe_allow_html=True)
        cols_per_row = 3
        for i in range(0, len(tools_list), cols_per_row):
            row_tools = tools_list[i:i+cols_per_row]; cols = st.columns(cols_per_row)
            for idx, tool in enumerate(row_tools):
                if idx < len(row_tools):
                    with cols[idx]:
                        with st.container():
                            st.markdown('<div class="tool-card-style">', unsafe_allow_html=True)
                            st.markdown(f'<h4 class="tool-card-title">{tool["name"]}</h4>', unsafe_allow_html=True); st.markdown(f'<p class="tool-card-description">{tool["description"]}</p>', unsafe_allow_html=True)
                            button_key = f"open_{tool['name'].replace(' ', '_').replace('/', '').replace('(', '').replace(')', '').lower()}_{title.replace(' ', '_').lower()}_{i+idx}"
                            st.markdown('<div class="tool-card-button-wrapper">', unsafe_allow_html=True)
                            if st.button(f"Open {tool['name']}", key=button_key, use_container_width=False):
                                if tool["iframe_able"]: st.session_state.metagen_tool_iframe_url = tool["url"]; st.session_state.metagen_tool_iframe_name = tool["name"]; st.rerun()
                                else: st.markdown(f"<script>window.open('{tool['url']}', '_blank');</script>", unsafe_allow_html=True); st.toast(f"Opening {tool['name']} in a new tab...", icon="🌐")
                            st.markdown('</div>', unsafe_allow_html=True)
                            st.markdown('</div>', unsafe_allow_html=True)
        st.markdown("<br>", unsafe_allow_html=True)
    if st.session_state.metagen_tool_iframe_url:
        st.subheader(f"Viewing: {st.session_state.metagen_tool_iframe_name}")
        if st.button(f"← Back to Metagenomics Resources", key="back_to_meta_resources_top_iframe"): st.session_state.metagen_tool_iframe_url = None; st.session_state.metagen_tool_iframe_name = None; st.rerun()
        try:
            st.components.v1.iframe(st.session_state.metagen_tool_iframe_url, height=800, scrolling=True)
            st.caption(f"Note: Viewing {st.session_state.metagen_tool_iframe_name} embedded. Functionality might be limited.")
        except Exception as e_iframe: st.error(f"Could not embed {st.session_state.metagen_tool_iframe_name}. Error: {e_iframe}"); st.markdown(f"Visit directly: <a href='{st.session_state.metagen_tool_iframe_url}' target='_blank'>{st.session_state.metagen_tool_iframe_url}</a>", unsafe_allow_html=True)
    else:
        for category_title, tools_in_category in tool_categories.items(): display_tool_section_cards(category_title, tools_in_category)
    st.markdown("---")

# --- MAIN APP STRUCTURE ---
PAGES = {
    "Home": render_home_page,
    "Genomic Toolkit": render_genomic_toolkit_page,
    "FASTA Analysis Utilities": render_fasta_analysis_tools_page,
    "Metagenomics Resources": render_metagenomics_tools_page,
    "About": render_about_section,
}

# Sidebar Navigation
st.sidebar.markdown("---")
page_selection = st.sidebar.radio("Go to", list(PAGES.keys()))
st.sidebar.markdown("<hr style='border-top: 1.5px solid #0a3a4a; margin: 1em 0;'>", unsafe_allow_html=True)


PAGES[page_selection]() 

# Sidebar footer
st.sidebar.markdown("---")
with st.sidebar.expander("GenVedix"):
    st.markdown(""" 
    **Last Updated-May 2025**  
    **Developed By Tejal Kale with ❤ for the Bioinformatics Community**  
    """)

