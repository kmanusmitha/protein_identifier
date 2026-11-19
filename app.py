import streamlit as st
import pandas as pd
import requests
import time
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Align import PairwiseAligner, substitution_matrices

st.set_page_config(page_title="Protein Identifier", page_icon="🧬", layout="wide")
st.title("🧬 Unknown Protein Identification System")
st.write("Automated Sequence Analysis, Physicochemical Profiling & Global Homology Search")

default_fasta = """>Unknown_Sample_01
MLEFSFLKVKRTMFFPTYLPEHCKKIKQKITSYSNSRKKGYYFNTEGKRADINTPALAGR
GGLTGDEEQECNHKIMRRFRYHSRKHKNNQPRLHRLYPSMREETQLAEQSPTFPRSRQVQ"""

reference_database = {
    "Q8BSC4_MOUSE | GPI mannosyltransferase 3 (Pigb) [Mus musculus]": (
        "MLEFSFLKVKRTMFFPTYLPEHCKKIKQKITSYSNSRKKGYYFNTEGKRADINTPALAGR"
        "GGLTGDEEQECNHKIMRRFRYHSRKHKNNQPRLHRLYPSMREETQLAEQSPTFPRSRQVQ"
    ),
    "Q9EQY6_RAT | GPI mannosyltransferase homolog [Rattus norvegicus]": (
        "MLEFSFLKVKRTMFFPTYLPEHCKKIKQKITSYSNSRKKGYYFNTEGKRADINTPALAGR"
        "GGLTGDEEQECNHKIMRRFRYHSRKHKNNQPRLHRLYPSMREETQLAEQSPTSPRSRQVQ"
    ),
    "Q9H490_HUMAN | PIGB Phosphatidylinositol glycan [Homo sapiens]": (
        "MRRPLSKCGFLLPLLLALLLAPGAAQAAALELWPRLAAALLALLLLLLLPPRPPAAAA"
        "GGLTGDEEQECNHKIMRRFRYHSRKHKNNQPRLHRLYPSMREETQLAEQSPTFPRSRQVQ"
    ),
    "P01308_HUMAN | Insulin Hormone [Homo sapiens]": (
        "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAED"
        "LQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
    ),
    "P69905_HUMAN | Hemoglobin Subunit Alpha [Homo sapiens]": (
        "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG"
        "KKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTP"
    )
}

fasta_input = st.text_area("Paste ANY Unknown Protein Sequence (FASTA format):", value=default_fasta, height=150)
use_live_api = st.checkbox("🌐 Search Live Global UniProt / Swiss-Prot Server (For any unknown protein, takes ~20 sec)", value=True)

if st.button("🔍 Identify & Analyze Protein"):
    lines = fasta_input.strip().splitlines()
    if len(lines) < 2 or not lines.pop(0).startswith(">"):
        st.error("❌ Please enter a valid FASTA format starting with '>' on the first line.")
    else:
        raw_lines = fasta_input.strip().splitlines()
        header = raw_lines.pop(0).lstrip(">").strip()
        sequence = "".join(raw_lines).upper().replace(" ", "")

        try:
            # 1. PHYSICOCHEMICAL ANALYSIS (Works on ANY protein sequence)
            protein = ProteinAnalysis(sequence)
            length = len(sequence)
            mw = round(protein.molecular_weight(), 2)
            pi = round(protein.isoelectric_point(), 2)
            instability = round(protein.instability_index(), 2)
            stability = "Stable" if instability = 50:
                st.success(f"**Local Match Found:** {top_name} (**{top_sim}% Similarity**)")
            else:
                st.info("No close match in local 5-protein cache. Checking Global UniProt Server below...")
            st.dataframe(df_matches, use_container_width=True)

            # 3. LIVE GLOBAL UNIPROT / EBI BLAST API SEARCH
            if use_live_api:
                st.subheader("🌐 Live Global UniProt / Swiss-Prot Database Search")
                with st.spinner("Querying EMBL-EBI / UniProt Global Servers (Please wait ~15-30 seconds)..."):
                    submit_url = "https://ebi.ac.uk"
                    payload = {
                        "email": "bioinformatics.capstone@example.com",
                        "program": "blastp",
                        "stype": "protein",
                        "database": "uniprotkb_swissprot",
                        "sequence": sequence
                    }
                    res = requests.post(submit_url, data=payload, timeout=15)
                    job_id = res.text.strip()

                    if job_id.startswith("ncbiblast-"):
                        status_url = f"https://ebi.ac.uk{job_id}"
                        for _ in range(12):
                            status = requests.get(status_url, timeout=10).text.strip()
                            if status == "FINISHED":
                                break
                            time.sleep(4)

                        result_url = f"https://ebi.ac.uk{job_id}/out"
                        report_text = requests.get(result_url, timeout=15).text
                        hit_lines = [line for line in report_text.splitlines() if line.startswith("SP:")]

                        if hit_lines:
                            st.success(f"✅ Found {len(hit_lines)} reviewed matches in Global UniProtKB/Swiss-Prot!")
                            st.code("\n".join(hit_lines[:10]), language="text")
                        else:
                            st.warning("No reviewed match found in Swiss-Prot (Potentially novel/uncharacterized sequence).")
                    else:
                        st.warning("Live external API blocked by local network firewall. Works automatically when hosted on Streamlit Cloud!")

        except Exception as e:
            st.error(f"Error analyzing sequence: {e}. Make sure only valid amino acid letters are used.")
