import sys
import os
import json
import base64
import streamlit as st
import matplotlib.pyplot as plt

# Download Graph imports
from io import BytesIO

def download_plot(fig, filename, key):
    buf = BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight")

    st.download_button(
        label="⬇️ Download Graph",
        data=buf.getvalue(),
        file_name=filename,
        mime="image/png",
        key=key
    )
    
# 🔬 Biopython imports
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO
import io

# Fix import path for ai_engine
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(BASE_DIR)

from ai_engine.analysis.gc_content import gc_content
from ai_engine.optimization.optimizer import optimize_codon
from ai_engine.optimization.harmonizer import harmonize_codon
from ai_engine.ml.predict import predict_expression
from ai_engine.utils.validator import detect_sequence_type
from ai_engine.analysis.gc_window import gc_content_window
from ai_engine.analysis.cai import calculate_cai
from ai_engine.analysis.codon_usage import codon_usage
from ai_engine.analysis.rare_codon import find_rare_codons
from ai_engine.analysis.ribosome import ribosome_profile
from ai_engine.analysis.mrna_structure import sliding_delta_g

# Constants
IMG_DIR = os.path.join(BASE_DIR, "images")
banner_path = os.path.join(IMG_DIR, "ai_codon.png")

# Streamlit Page config
st.set_page_config(
    page_title="AI Codon Optimizer",
    page_icon="🧬",
    layout="wide"
)

# ====== Style: Fresh look & page-specific backgrounds ======
st.markdown("""
<style>
/* Full-page gradient background */
body {
    background: linear-gradient(135deg, #f0f4f8 0%, #e2ebf0 100%);
    color: #111827;
}

/* Hero banner */
.hero {
    background-size: cover;
    padding: 60px;
    border-radius: 15px;
    color: white;
}

/* Card styling */
div[class*="stCard"] {
    background-color: rgba(255, 255, 255, 0.95);
    border-radius: 15px;
    box-shadow: 0 4px 20px rgba(0,0,0,0.1);
    padding: 20px;
}

/* Streamlit buttons */
.stButton>button {
    background-color: #1e3a8a;
    color: white;
    border-radius: 8px;
    padding: 8px 20px;
    font-weight: bold;
}

/* Text areas */
textarea {
    border-radius: 8px;
    border: 1px solid #ccc;
    padding: 8px;
    font-family: monospace;
}

/* Sidebar gradient */
section[data-testid="stSidebar"] {
    background: linear-gradient(180deg, #0b1b3a 0%, #02112e 100%);
}

/* Sidebar text */
section[data-testid="stSidebar"] * {
    color: white !important;
}
</style>
""", unsafe_allow_html=True)

# ====== Banner Image ======
def get_base64(img_path):
    if os.path.exists(img_path):
        with open(img_path, "rb") as f:
            return base64.b64encode(f.read()).decode()
    return ""

banner_base64 = get_base64(banner_path).replace("\n", "")

# Sidebar navigation styling and menu
with st.sidebar:
    st.markdown(
        """
        <style>
        [data-testid="stSidebar"] {
            padding-top: 20px;
        }
        .sidebar-text {
            color: white !important;
            font-size: 20px;
            font-weight: 600;
            padding: 10px 12px;
            display: flex;
            align-items: center;
            gap: 8px;
            white-space: nowrap;
        }
        .stRadio > div {
            padding-left: 10px;
        }
        .stRadio > div > label {
            cursor: pointer;
        }
        </style>
        """,
        unsafe_allow_html=True
    )
    menu = st.radio(
        "",
        ["🏠 Home", "📘 About", "📞 Contact"],
        label_visibility="collapsed",
        index=0,
        key="menu",
    )

# Load codon tables
codon_dir = os.path.join(BASE_DIR, "ai_engine", "codon_tables")
if not os.path.exists(codon_dir):
    os.makedirs(codon_dir)

organisms = [
    f.replace(".json", "")
    for f in os.listdir(codon_dir)
    if f.endswith(".json")
]

# ======== PAGE CONTENT ========

if menu == "🏠 Home":
    # Page-specific background
    st.markdown("""
    <style>body {background: linear-gradient(to right, #f8fafc, #e0f2fe);}</style>
    """, unsafe_allow_html=True)

    # HERO
    st.markdown(f"""
    <div class="hero" style="background-image: linear-gradient(rgba(15,23,42,0.75), rgba(30,58,138,0.75)),
                      url('data:image/jpeg;base64,{banner_base64}')">
        <h1>🧬 AI Codon Optimizer</h1>
        <p>Optimize gene sequences using AI & Bioinformatics</p>
    </div>
    """, unsafe_allow_html=True)
    
    st.markdown('<div class="card">', unsafe_allow_html=True)

    # FASTA Upload
    uploaded_fasta = st.file_uploader("Upload FASTA File", type=["fasta", "fa", "txt"])
    sequence = ""

    def format_sequence(seq, line_length=60):
        return "\n".join([seq[i:i+line_length] for i in range(0, len(seq), line_length)])

    if uploaded_fasta:
        try:
            fasta_data = uploaded_fasta.read().decode("utf-8")
            records = list(SeqIO.parse(StringIO(fasta_data), "fasta"))
            if len(records) == 0:
                st.error("❌ No valid sequences found in FASTA")
            else:
                if len(records) > 1:
                    st.warning(f"⚠️ Multiple sequences found. Using first: {records[0].id}")
                sequence = str(records[0].seq)
                st.success(f"✅ Loaded: {records[0].id}")
                formatted_seq = format_sequence(sequence)
                st.subheader("Sequence")
                st.text_area("", formatted_seq, height=100)
                st.caption(f"Length: {len(sequence)} bases")
                st.download_button("Download Sequence", data=sequence, file_name="sequence.txt")
        except Exception as e:
            st.error(f"Error reading FASTA: {e}")

    # Manual Input
    manual_seq = st.text_area("Enter DNA / Protein Sequence")
    if manual_seq:
        sequence = manual_seq.strip()

    # Organism
    organism = st.selectbox("Select Organism", organisms) if organisms else None

    # Mode
    mode = st.radio(
        "Select Optimization Mode",
        ["Max Expression (Optimization)", "Balanced Folding (Harmonization)"]
    )

    # ANALYZE
    if st.button("Analyze"):
        if sequence:
            seq_type = detect_sequence_type(sequence)
            st.write("Sequence Type:", seq_type)
            if seq_type != "INVALID":
                dna_seq = sequence if seq_type == "DNA" else optimize_codon(sequence, organism)
                if dna_seq:
                    gc = gc_content(dna_seq)
                    st.write("GC Content:", gc)
                    gc_vals = gc_content_window(dna_seq)
                    fig, ax = plt.subplots(figsize=(6, 2.5))
                    ax.plot(gc_vals)
                    ax.set_title("GC Content (Sliding Window)")
                    ax.set_xlabel("Position")
                    ax.set_ylabel("GC Ratio")
                    st.pyplot(fig, use_container_width=False)
                    download_plot(fig, "gc_content.png", key="gc_download")
                    plt.close(fig)

    # OPTIMIZE
    if st.button("Optimize"):
        if sequence and organism:
            original_seq = sequence
            optimized = optimize_codon(sequence, organism) if mode.startswith("Max") else harmonize_codon(sequence, organism)
            if not optimized:
                st.error("Optimization failed")
                st.stop()
            st.write("🧬 Result:", optimized)
            # Sequence Comparison
            st.subheader("🔍 Sequence Comparison")
            col1, col2 = st.columns(2)
            with col1:
                st.markdown("**Original Sequence**")
                st.text_area("", format_sequence(original_seq), height=150)
            with col2:
                st.markdown("**Optimized Sequence**")
                st.text_area("", format_sequence(optimized), height=150)
            cai_score = calculate_cai(optimized)
            st.write("🧠 CAI Score:", cai_score)
            # Codon Usage
            usage = codon_usage(optimized)
            fig, ax = plt.subplots(figsize=(7, 3))
            ax.bar(list(usage.keys()), list(usage.values()))
            plt.xticks(rotation=90)
            ax.set_title("Codon Usage")
            st.pyplot(fig, use_container_width=False)
            download_plot(fig, "codon_usage.png", key="codon_download")
            plt.close(fig)
            # Rare Codons
            rare = find_rare_codons(optimized)
            if rare:
                st.warning(f"{len(rare)} rare codons found")
            # Ribosome Profile
            speeds = ribosome_profile(optimized)
            fig3, ax3 = plt.subplots(figsize=(6, 2.5))
            ax3.plot(speeds)
            ax3.set_title("Ribosome Speed")
            st.pyplot(fig3, use_container_width=False)
            download_plot(fig3, "ribosome_profile.png", key="ribosome_download")
            plt.close(fig3)
            # mRNA ΔG
            dg_vals = sliding_delta_g(optimized)
            fig4, ax4 = plt.subplots(figsize=(6, 2.5))
            ax4.plot(dg_vals)
            ax4.set_title("mRNA ΔG")
            st.pyplot(fig4, use_container_width=False)
            download_plot(fig4, "mrna_dg.png", key="dg_download")
            plt.close(fig4)
            # Warnings
            if cai_score < 0.7:
                st.warning("⚠️ Low expression potential → consider optimization")
            if len(rare) > 5:
                st.warning("⚠️ Many rare codons detected → may slow translation")
            if dg_vals and min(dg_vals) < -30:
                st.warning("⚠️ Strong mRNA structure → may block ribosome")
            # Translation
            protein = Seq(optimized).translate()
            st.write("Protein:", protein)
            # Download FASTA
            record = SeqRecord(Seq(optimized), id="Optimized")
            fasta_io = io.StringIO()
            SeqIO.write(record, fasta_io, "fasta")
            st.download_button("Download FASTA", data=fasta_io.getvalue(), file_name="optimized.fasta")

    # AI Prediction
    if st.button("Predict Expression Score"):
        if sequence:
            st.write("Score:", predict_expression(sequence))

    # Upload Custom Table
    uploaded_file = st.file_uploader("Upload custom codon table (JSON)")
    if uploaded_file:
        custom_table = json.load(uploaded_file)
        with open(os.path.join(codon_dir, "custom.json"), "w") as f:
            json.dump(custom_table, f)
        st.success("Custom organism added!")

elif menu == "📘 About":
    st.markdown("""<style>body {background-color: #fef9c3;}</style>""", unsafe_allow_html=True)
    st.title("📘 About")
    st.write("""
The AI Codon Optimizer is an advanced bioinformatics web application developed to enhance gene sequences for improved protein expression across different host organisms. It combines classical codon optimization strategies with modern computational and AI-driven approaches to provide a comprehensive sequence engineering platform. The application accepts both DNA and protein sequences as input and intelligently processes them by identifying sequence type, converting protein sequences into optimized DNA, and applying organism-specific codon usage patterns. Users can choose between maximum expression optimization, which prioritizes highly preferred codons, and codon harmonization, which preserves translational kinetics for proper protein folding.

Beyond optimization, the platform performs a wide range of analytical evaluations to assess sequence quality and expression potential. These include:

. GC Content Analysis
. Codon Adaptation Index (CAI) for expression efficiency
. Codon Usage Profiling
. Rare Codon Identification
. Ribosome Translation Speed Profiling
. mRNA Secondary Structure Prediction (ΔG analysis)

To further enhance usability, the tool integrates a machine learning–based module to predict gene expression levels, providing users with an additional layer of insight. It also generates suggestions and warnings, helping users identify potential bottlenecks such as strong secondary structures, rare codon clusters, or suboptimal GC content. This platform supports custom codon tables, allowing flexibility for non-model organisms and specialized research needs. Additionally, users can export optimized sequences and generate comprehensive reports, making it suitable for documentation, collaboration, and downstream experimental workflows.

Overall, the AI Codon Optimizer serves as an integrated solution for researchers, students, and professionals working in Biotechnology, Synthetic biology, Pharmaceutical research, Genetic engineering. By combining automation, visualization, and AI-driven insights, the tool simplifies complex sequence optimization tasks and enables more informed decision-making in gene design and expression studies.
""")

elif menu == "📞 Contact":
    st.markdown("""<style>body {background-color: #fde2e4;}</style>""", unsafe_allow_html=True)
    st.title("📞 Contact Us")
    col1, col2 = st.columns(2)
    st.markdown("""
    <style>
    .profile-card { text-align:center; }
    .profile-card img {
        width:220px;
        height:220px;
        object-fit:cover;
        border-radius:15px;
        box-shadow:0 4px 15px rgba(0,0,0,0.15);
    }
    </style>
    """, unsafe_allow_html=True)

    def img64(path):
        return base64.b64encode(open(path, "rb").read()).decode()

    tejas = os.path.join(IMG_DIR, "tejas.jpeg")
    mohasin = os.path.join(IMG_DIR, "Mohasin.jpeg")

    with col1:
        if os.path.exists(tejas):
            st.markdown(f"""
            <div class="profile-card">
                <img src="data:image/jpeg;base64,{img64(tejas)}">
                <h4>Tejasman Kour</h4>
                <p>Bioinformatics Developer</p>
                <p>Linkedin: linkedin.com/in/tejasman-kour-705659255</p>
                <p>📧 tejasmankouranand@gmail.com</p>
            </div>
            """, unsafe_allow_html=True)

    with col2:
        if os.path.exists(mohasin):
            st.markdown(f"""
            <div class="profile-card">
                <img src="data:image/jpeg;base64,{img64(mohasin)}">
                <h4>Mohasin Sayyad</h4>
                <p>Project Engineer & UI Developer</p>
                <p>Linkedin: linkedin.com/in/mohasinsayyad31294</p>
                <p>📧 mohasins@cdac.in</p>
            </div>
            """, unsafe_allow_html=True)

    st.markdown("---")
    st.write("📧 support.codon.ai@gmail.com")
