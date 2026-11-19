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
    if not lines or len(lines) == 1 or not lines.pop(0).startswith(">"):
        st.error("❌ Please enter a valid FASTA format starting with '>' on the first line.")
    else:
        raw_lines = fasta_input.strip().splitlines()
        header = raw_lines.pop(0).lstrip(">").strip()
        sequence = "".join(raw_lines).upper().replace(" ", "")

        try:
            protein = ProteinAnalysis(sequence)
            length = len(sequence)
            mw = round(protein.molecular_weight(), 2)
            pi = round(protein.isoelectric_point(), 2)
            instability = round(protein.instability_index(), 2)
            if int(instability) in range(0, 40):
                stability = "Stable"
            else:
                stability = "Unstable"

            st.subheader(f"📊 Analysis Report for: `{header}`")
            col1, col2, col3, col4 = st.columns(4)
            col1.metric("Length", f"{length} A.A.")
            col2.metric("Molecular Weight", f"{mw} Da")
            col3.metric("Isoelectric Point (pI)", pi)
            col4.metric("Instability Index", f"{instability} ({stability})")

            st.subheader("🧪 Amino Acid Composition")
            aa_counts = protein.count_amino_acids()
            df_aa = pd.DataFrame(list(aa_counts.items()), columns=["Amino Acid", "Count"]).set_index("Amino Acid")
            st.bar_chart(df_aa)

            st.subheader("🧬 Local Reference Homology (BLOSUM62)")
            aligner = PairwiseAligner()
            aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
            aligner.open_gap_score = -10
            aligner.extend_gap_score = -0.5

            matches_list = []
            for name, ref_seq in reference_database.items():
                score = aligner.score(sequence, ref_seq)
                best_alignment = next(iter(aligner.align(sequence, ref_seq)))
                identical = best_alignment.counts().identities
                identity_pct = round((identical / len(sequence)) * 100, 2)
                matches_list.append({
                    "Matched Reference Protein": name,
                    "Similarity (%)": identity_pct,
                    "BLOSUM62 Score": round(score, 1)
                })

            df_matches = pd.DataFrame(matches_list).sort_values(by="BLOSUM62 Score", ascending=False)
            top_name = df_matches["Matched Reference Protein"].tolist().pop(0)
            top_sim = df_matches["Similarity (%)"].tolist().pop(0)

            if top_sim == 0:
                st.info("No close match in local 5-protein cache. Checking Global UniProt Server below...")
            else:
                st.success(f"**Local Match Found:** {top_name} (**{top_sim}% Similarity**)")
            st.dataframe(df_matches, use_container_width=True)

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
                        st.warning("Live external API busy. Please try again in a moment.")

        except Exception as e:
            st.error(f"Error analyzing sequence: {e}. Make sure only valid amino acid letters are used.")
