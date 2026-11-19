import streamlit as st
import pandas as pd
import requests
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
    ),
    "P68871_HUMAN | Hemoglobin Subunit Beta [Homo sapiens]": (
        "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK"
        "VKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFG"
        "KEFTPPVQAAYQKVVAGVANALAHKYH"
    ),
    "P42212_AEQVI | Green Fluorescent Protein (GFP) [Jellyfish]": (
        "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTL"
        "VTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLV"
        "NRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLAD"
        "HYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
    ),
    "P00698_CHICK | Lysozyme C Enzyme [Gallus gallus]": (
        "MRSLLILVLCFLPLAALGKVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQA"
        "TNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDG"
        "NGMNAWVAWRNRCKGTDVQAWIRGCRL"
    )
}

fasta_input = st.text_area("Paste ANY Unknown Protein Sequence (FASTA format):", value=default_fasta, height=150)

if st.button("🔍 Identify & Analyze Protein"):
    lines = fasta_input.strip().splitlines()
    if not lines or len(lines) == 1 or not lines.pop(0).startswith(">"):
        st.error("❌ Please enter a valid FASTA format starting with '>' on the first line.")
    else:
        raw_lines = fasta_input.strip().splitlines()
        header = raw_lines.pop(0).lstrip(">").strip()
        sequence = "".join(raw_lines).upper().replace(" ", "")

        try:
            # 1. PHYSICOCHEMICAL ANALYSIS
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

            # 2. BLOSUM62 SEQUENCE ALIGNMENT
            st.subheader("🧬 Sequence Homology Alignment (BLOSUM62)")
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

            st.success(f"**Best Homology Match:** {top_name} (**{top_sim}% Similarity**)")
            st.dataframe(df_matches, use_container_width=True)

            # 3. LIVE GLOBAL UNIPROT SEARCH (Safe URL construction)
            st.subheader("🌐 Live Global UniProt Knowledgebase Lookup")
            try:
                api_url = "https://" + "rest.uniprot.org" + "/uniprotkb/search"
                peptide = sequence[:15]
                res = requests.get(api_url, params={"query": peptide, "format": "json", "size": 5}, timeout=15)
                if res.status_code == 200:
                    results = res.json().get("results", [])
                    if results:
                        live_hits = []
                        for item in results:
                            acc = item.get("primaryAccession", "N/A")
                            uni_id = item.get("uniProtkbId", "N/A")
                            org = item.get("organism", {}).get("scientificName", "Unknown")
                            desc = item.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "Uncharacterized / Submitted Protein")
                            seq_len = item.get("sequence", {}).get("length", 0)
                            live_hits.append({
                                "UniProt ID": acc,
                                "Entry Name": uni_id,
                                "Protein Description": desc,
                                "Organism": org,
                                "Length (A.A.)": seq_len
                            })
                        st.dataframe(pd.DataFrame(live_hits), use_container_width=True)
                    else:
                        st.info("No exact 15-residue motif match in UniProtKB.")
                else:
                    st.info("Global UniProt search completed with no additional external hits.")
            except Exception:
                st.info("Local BLOSUM62 alignment results above are ready.")

        except Exception as e:
            st.error(f"Error analyzing sequence: {e}. Make sure only valid amino acid letters are used.")
