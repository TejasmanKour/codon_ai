# 🧬 AI Codon Optimizer

A powerful **bioinformatics + AI-based web application** that optimizes gene sequences to improve protein expression in target organisms.

---

## 🚀 Features

* 🔍 **Automatic Sequence Detection** (DNA / Protein)
* 🧬 **Codon Optimization Engine**
* ⚖️ **Harmonization Mode for Balanced Folding**
* 📊 **GC Content Analysis**
* 🧠 **Codon Adaptation Index (CAI) Calculation**
* 📈 **Codon Usage Visualization**
* ⚠️ **Rare Codon Detection**
* 🔬 **mRNA Structure Analysis (ΔG)**
* 🤖 **AI-Based Expression Prediction**
* 🔍 **Side-by-Side Sequence Comparison**
* 📥 **Download Graphs & FASTA Output**
* 🧑‍💻 **Multi-page UI (Home, About, Contact)**

---

## 🧪 How It Works

1. Upload a FASTA file or enter a DNA/Protein sequence
2. System detects sequence type automatically
3. Converts protein → DNA (if required)
4. Applies organism-specific codon optimization
5. Performs biological analysis (GC, CAI, etc.)
6. Displays graphs, warnings, and optimized sequence
7. Allows download of results

---

## 🛠️ Tech Stack

* **Frontend/UI:** Streamlit
* **Backend:** Python
* **Bioinformatics:** Biopython
* **Visualization:** Matplotlib
* **AI/ML:** Scikit-learn (Expression Prediction)

---

## 📁 Project Structure

```
codon_ai_streamlit/
│
├── app.py
├── ai_engine/
│   ├── analysis/
│   ├── optimization/
│   ├── utils/
│   ├── ml/
│   └── codon_tables/
├── images/
└── requirements.txt
```

---

## ⚙️ Installation & Setup

### 1️⃣ Clone or Copy Project

```
git clone <your-repo-url>
cd codon_ai_streamlit
```

### 2️⃣ Create Virtual Environment

```
python -m venv venv
```

Activate:

* Windows:

```
venv\Scripts\activate
```

* Linux/Mac:

```
source venv/bin/activate
```

---

### 3️⃣ Install Dependencies

```
pip install -r requirements.txt
```

(If not available:)

```
pip install streamlit biopython matplotlib numpy scikit-learn
```

---

### 4️⃣ Run the Application

```
streamlit run app.py
```

Open in browser:

```
http://localhost:8501
```

---

## 📊 Output

* Optimized DNA Sequence
* Protein Translation
* Graphs (GC Content, Codon Usage, etc.)
* Expression Score Prediction
* Downloadable FASTA file

---

## 📸 Screenshots

> Add screenshots of your app here (Home, Graphs, Comparison, Contact page)

---

## 🎯 Applications

* 🧪 Biotechnology Research
* 💊 Pharmaceutical Industry
* 🧬 Synthetic Biology
* 🦠 Microbial Engineering
* 🎓 Academic Learning
* 🤖 AI in Bioinformatics

---

## 👥 Team

* **Mohasin Sayyad** – Bioinformatics Developer
* **Teammate 2** – AI/ML Engineer
* **Teammate 3** – Backend Developer

---

## 📞 Contact

📧 Email: [support@codonai.com](mailto:support@codonai.com)
📍 Location: India

---

## ⭐ Future Improvements

* 🔥 Mutation highlighting (color diff)
* 📊 Expression improvement %
* 🎨 Dark mode toggle
* 🧬 Advanced ML models

---

## 📜 License

This project is developed for academic purposes.
Feel free to use and modify.

---

## 🙌 Acknowledgement

This project combines **bioinformatics, AI, and web technologies** to solve real-world biological problems in gene optimization.
