#!/usr/bin/env python3
"""
SMILES to Oligo Notation — Streamlit Web App
=============================================
Wraps smiles_to_oligo_notation3.py into a clean web interface.

Run:  streamlit run app.py
"""

import streamlit as st
from smiles_to_oligo_notation3 import smiles_to_all

# ============================================================
# Page config
# ============================================================
st.set_page_config(
    page_title="SMILES → Oligo Notation",
    page_icon="🧬",
    layout="wide",
)

# ============================================================
# Custom CSS — Apple-designer pastel style
# ============================================================
st.markdown("""
<style>
    /* ---- Global ---- */
    .stApp {
        background: linear-gradient(135deg, #f5f0ff 0%, #eaf6ff 50%, #f0fff4 100%);
    }
    section[data-testid="stSidebar"] {
        background: #f8f5ff;
    }

    /* ---- Title area ---- */
    .main-title {
        text-align: center;
        font-size: 2.6rem;
        font-weight: 700;
        color: #4a3580;
        margin-bottom: 0;
        letter-spacing: -0.5px;
    }
    .sub-title {
        text-align: center;
        font-size: 1.05rem;
        color: #7c6ea6;
        margin-top: 0;
        margin-bottom: 2rem;
    }

    /* ---- Cards ---- */
    .result-card {
        background: white;
        border-radius: 16px;
        padding: 1.5rem 1.8rem;
        margin-bottom: 1.2rem;
        box-shadow: 0 2px 12px rgba(120, 100, 180, 0.08);
        border: 1px solid #ede8f5;
    }
    .card-header {
        font-size: 1.1rem;
        font-weight: 600;
        color: #5b4a8a;
        margin-bottom: 0.6rem;
        display: flex;
        align-items: center;
        gap: 8px;
    }
    .card-label {
        font-size: 0.85rem;
        font-weight: 500;
        color: #9080b8;
        margin-bottom: 0.25rem;
    }
    .seq-box {
        background: #faf8ff;
        border: 1px solid #e6e0f3;
        border-radius: 10px;
        padding: 0.8rem 1rem;
        font-family: 'SF Mono', 'Fira Code', 'Consolas', monospace;
        font-size: 0.88rem;
        word-wrap: break-word;
        white-space: pre-wrap;
        color: #3d2e6b;
        margin-bottom: 0.8rem;
        line-height: 1.6;
    }

    /* ---- MW table ---- */
    .mw-table {
        width: 100%;
        border-collapse: separate;
        border-spacing: 0;
        border-radius: 10px;
        overflow: hidden;
        margin-top: 0.5rem;
    }
    .mw-table th {
        background: #ede8f5;
        color: #5b4a8a;
        text-align: left;
        padding: 0.6rem 1rem;
        font-size: 0.85rem;
        font-weight: 600;
    }
    .mw-table td {
        padding: 0.55rem 1rem;
        font-size: 0.88rem;
        color: #4a4070;
        border-bottom: 1px solid #f0ecf7;
    }
    .mw-table tr:last-child td {
        border-bottom: none;
    }

    /* ---- Conjugate badge ---- */
    .conj-badge {
        display: inline-block;
        background: linear-gradient(135deg, #fff0f5, #f5e6ff);
        border: 1px solid #e6d0f0;
        border-radius: 8px;
        padding: 0.5rem 0.9rem;
        margin-bottom: 0.5rem;
        font-size: 0.88rem;
        color: #6b4580;
    }

    /* ---- Misc ---- */
    .stTextArea > div > div > textarea {
        border-radius: 12px !important;
        border: 1.5px solid #d8d0e8 !important;
        font-family: 'SF Mono', 'Fira Code', monospace !important;
        font-size: 0.85rem !important;
    }
    div.stButton > button {
        background: linear-gradient(135deg, #7c5cbf, #9b7ed8);
        color: white;
        border: none;
        border-radius: 10px;
        padding: 0.55rem 2.5rem;
        font-size: 1rem;
        font-weight: 600;
        letter-spacing: 0.3px;
        transition: all 0.2s;
    }
    div.stButton > button:hover {
        background: linear-gradient(135deg, #6a4aad, #8a6cc8);
        box-shadow: 0 4px 15px rgba(124, 92, 191, 0.3);
    }
    .success-banner {
        background: linear-gradient(135deg, #e8ffe8, #f0fff4);
        border: 1px solid #c3e6c3;
        border-radius: 10px;
        padding: 0.5rem 1rem;
        color: #2e7d32;
        font-size: 0.9rem;
        margin-bottom: 1rem;
        text-align: center;
    }
    .footer-text {
        text-align: center;
        font-size: 0.78rem;
        color: #b0a0c8;
        margin-top: 3rem;
    }
</style>
""", unsafe_allow_html=True)


# ============================================================
# Header
# ============================================================
st.markdown('<p class="main-title">🧬 SMILES → Oligo Notation</p>',
            unsafe_allow_html=True)
st.markdown(
    '<p class="sub-title">Paste a SMILES string to get Thermo &amp; Agilent '
    'sequences, molecular weights, and conjugate analysis</p>',
    unsafe_allow_html=True,
)

# ============================================================
# Input area
# ============================================================
col_input, col_info = st.columns([3, 1])

with col_input:
    smiles_input = st.text_area(
        "Enter SMILES",
        height=120,
        placeholder="Paste SMILES here (single strand or duplex with '.' separator)…",
        label_visibility="collapsed",
    )

with col_info:
    st.markdown("""
    <div class="result-card" style="padding:1rem 1.2rem;">
        <div class="card-header">💡 Tips</div>
        <ul style="font-size:0.82rem; color:#6b5a90; padding-left:1.2rem; margin:0;">
            <li>Single strand: paste one SMILES</li>
            <li>Duplex: separate strands with <code>.</code></li>
            <li>Supports PS, PO, VP, MOE, ADS, and more</li>
        </ul>
    </div>
    """, unsafe_allow_html=True)

# Analyze button
_, btn_col, _ = st.columns([2, 1, 2])
with btn_col:
    analyze = st.button("🔬  Analyze", use_container_width=True)

# ============================================================
# Run analysis
# ============================================================
if analyze and smiles_input and smiles_input.strip():
    try:
        thermo, agilent, mw_lines, conj_lines, alert_lines = smiles_to_all(smiles_input.strip())

        st.markdown('<div class="success-banner">✅ Analysis complete!</div>',
                    unsafe_allow_html=True)

        # ---- Sequences ----
        left_col, right_col = st.columns(2)

        with left_col:
            st.markdown(f"""
            <div class="result-card">
                <div class="card-header">🔴 Thermo BioPharma Finder</div>
                <div class="card-label">Full Thermo notation (base + sugar + linker)</div>
                <div class="seq-box">{thermo}</div>
            </div>
            """, unsafe_allow_html=True)

        with right_col:
            st.markdown(f"""
            <div class="result-card">
                <div class="card-header">🔵 Agilent BioConfirm</div>
                <div class="card-label">Standard Agilent notation</div>
                <div class="seq-box">{agilent}</div>
            </div>
            """, unsafe_allow_html=True)

        # ---- Copy-friendly code blocks ----
        with st.expander("📋 Copy sequences (plain text)"):
            st.code(f"Thermo:  {thermo}", language=None)
            st.code(f"Agilent: {agilent}", language=None)

        # ---- Molecular Weight ----
        mw_html = '<div class="result-card">'
        mw_html += '<div class="card-header">⚖️ Molecular Weight Summary</div>'
        mw_html += '<table class="mw-table">'
        mw_html += '<tr><th>Component</th><th>Monoisotopic (Da)</th>'
        mw_html += '<th>Average (Da)</th></tr>'

        current_label = ""
        mono_val = ""
        avg_val = ""
        for line in mw_lines:
            line = line.strip()
            if not line:
                continue
            if line.startswith("Monoisotopic:"):
                mono_val = line.split(":", 1)[1].strip()
            elif line.startswith("Average:"):
                avg_val = line.split(":", 1)[1].strip()
                mw_html += (f'<tr><td>{current_label}</td>'
                            f'<td>{mono_val}</td><td>{avg_val}</td></tr>')
                mono_val = ""
                avg_val = ""
            else:
                current_label = line

        mw_html += '</table></div>'
        st.markdown(mw_html, unsafe_allow_html=True)

        # ---- Conjugates ----
        conj_html = '<div class="result-card">'
        conj_html += '<div class="card-header">🔗 Conjugate Detection</div>'

        if conj_lines and conj_lines[0].strip() == "None detected":
            conj_html += ('<div class="conj-badge">'
                          '✅ No conjugates detected — standard oligonucleotide'
                          '</div>')
        else:
            for line in conj_lines:
                line = line.strip()
                if not line:
                    continue
                if line.startswith("Strand") or line.startswith("Position"):
                    conj_html += (f'<div class="conj-badge" '
                                  f'style="margin-top:0.5rem;">'
                                  f'<b>{line}</b></div>')
                elif line.startswith("Description:"):
                    conj_html += (f'<div style="margin-left:1rem; '
                                  f'font-size:0.85rem; color:#7a6090;">'
                                  f'{line}</div>')
                elif line.startswith("Monoisotopic:"):
                    conj_html += (f'<div style="margin-left:1rem; '
                                  f'font-size:0.85rem; color:#5a4a70;">'
                                  f'{line}</div>')
                elif line.startswith("Average:"):
                    conj_html += (f'<div style="margin-left:1rem; '
                                  f'font-size:0.85rem; color:#5a4a70; '
                                  f'margin-bottom:0.5rem;">{line}</div>')

        conj_html += '</div>'
        st.markdown(conj_html, unsafe_allow_html=True)
# ---- Modification Check ----
        alert_html = '<div class="result-card">'
        alert_html += '<div class="card-header">🔍 Modification Check</div>'
        for line in alert_lines:
            ls = line.strip()
            if not ls:
                continue
            if 'All modifications recognized' in ls:
                alert_html += '<div class="conj-badge" style="background:linear-gradient(135deg,#e8ffe8,#f0fff4);border-color:#c3e6c3;color:#2e7d32;">' + ls + '</div>'
            elif 'New/Unknown' in ls:
                alert_html += '<div style="color:#d32f2f;font-weight:600;font-size:0.95rem;">' + ls + '</div>'
            elif ls.startswith('['):
                alert_html += '<div style="background:#fef2f2;border:1px solid #fecaca;border-radius:8px;padding:0.5rem 0.8rem;margin:0.4rem 0;font-size:0.85rem;color:#7f1d1d;">' + ls + '</div>'
            elif ls.startswith('Symbol:'):
                alert_html += '<div style="margin-left:1rem;font-size:0.82rem;color:#991b1b;font-weight:500;">' + ls + '</div>'
            elif ls.startswith('->'):
                alert_html += '<div style="margin-left:1rem;font-size:0.82rem;color:#7a6090;">' + ls + '</div>'
            elif 'Tip' in ls:
                alert_html += '<div style="font-style:italic;color:#9080b8;font-size:0.82rem;margin-top:0.5rem;">' + ls + '</div>'
        alert_html += '</div>'
        st.markdown(alert_html, unsafe_allow_html=True)        

    except Exception as e:
        st.error(f"❌ Error analyzing SMILES: {str(e)}")

elif analyze:
    st.warning("⚠️ Please paste a SMILES string first.")

# ============================================================
# Footer
# ============================================================
st.markdown(
    '<p class="footer-text">Built with ❤️ by April Xia · '
    'LGM Analytical Team · Powered by RDKit + Streamlit</p>',
    unsafe_allow_html=True,
)