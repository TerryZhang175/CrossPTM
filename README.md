# CrossPTM: Structural PTM Crosstalk Analysis

**CrossPTM** is a bioinformatics pipeline and web interface designed to analyze the spatial crosstalk between a primary Post-Translational Modification (PTM) (e.g., Glutathionylation) and various other PTM types (e.g., Phosphorylation, Acetylation).

By leveraging **AlphaFold** predicted structures, CrossPTM calculates 3D Euclidean distances between modified residues to identify potential functional interplay and statistical enrichment.

## ⚠️Attention 
The Database/ folder contains only example files illustrating the required data format.
Please download the appropriate PTM database from dbPTM before running the program.

## 🚀 Features

*   **Interactive Web UI**: A modern, user-friendly interface to manage the entire analysis workflow.
*   **3D Structural Analysis**: Automatically fetches AlphaFold PDB structures and computes spatial distances between PTM sites.
*   **Parallelized Data Fetching**: Multi-threaded downloader for AlphaFold structures (8x concurrency) to speed up setup.
*   **Enrichment Analysis**: Statistical evaluation of PTM colocalization (Hypergeometric tests, Odds Ratios).
*   **Customizable**: 
    *   Upload any primary PTM dataset.
    *   Support for custom secondary PTM datasets.
    *   Configurable residue filtering (e.g., target only Cysteines or S/T/Y).
*   **Visual Reports**: Generates distance distribution summaries and enrichment reports.

## 🛠️ Installation

### Prerequisites
*   Python 3.9+
*   Internet connection (for fetching AlphaFold structures)

### Setup
1.  Clone the repository:
    ```bash
    git clone https://github.com/yourusername/CrossPTM.git
    cd CrossPTM
    ```

2.  Create a virtual environment (recommended):
    ```bash
    python3 -m venv .venv
    source .venv/bin/activate
    ```

3.  Install dependencies:
    ```bash
    pip install -r requirements.txt
    ```

## 🖥️ Usage

### Method 1: Web Interface (Recommended)
The web UI is the easiest way to run the pipeline.

1.  **Start the server**:
    ```bash
    # Local access only (default)
    uvicorn app.main:app --reload

    # Allow LAN access (accessible by others on your network)
    uvicorn app.main:app --host 0.0.0.0 --port 8000 --reload
    ```
2.  **Open your browser**:
    *   **Local**: [http://127.0.0.1:8000](http://127.0.0.1:8000)
    *   **LAN**: `http://<YOUR_IP_ADDRESS>:8000` (e.g., `http://192.168.1.5:8000`)

3.  **Follow the steps**:
    *   **Step 1**: Upload your Primary PTM file (CSV). Specify columns for UniProt ID and Site position.
    *   **Step 2**: The system will automatically download missing AlphaFold structures in the background.
    *   **Step 3**: Select existing PTM datasets (Phosphorylation, Ubiquitination, etc.) or upload your own Custom PTMs.
    *   **Step 4**: Click **Run Pipeline**.
    *   **Step 5**: View results, including pair counts, distance stats, and enrichment summaries.

### Method 2: Command Line Interface (Advanced)
You can also run individual scripts manually if you prefer the CLI.

1.  **Load Data**:
    ```bash
    python scripts/ptm_loader.py --glut-file path/to/primary.csv --all-ptms
    ```
2.  **Compute Distances**:
    ```bash
    python scripts/compute_distances.py --glut path/to/primary.csv
    ```
3.  **Run Enrichment**:
    ```bash
    python scripts/enrichment.py --glut path/to/primary.csv
    ```

## 📂 Project Structure

```text
CrossPTM/
├── app/                        # Web Application (FastAPI)
│   ├── main.py                 # API Entry point
│   ├── pipeline_runner.py      # Logic linking UI to scripts
│   └── static/                 # Frontend assets (HTML/CSS)
├── scripts/                    # Core Analysis Logic
│   ├── fetch_alphafold_structures.py  # Multi-threaded PDB downloader
│   ├── ptm_loader.py           # Data parsing and normalization
│   ├── compute_distances.py    # 3D distance calculations
│   └── enrichment.py           # Statistical analysis
├── Dataset/                    # Built-in PTM reference datasets
├── alphafold_structures/       # Downloaded PDB files (Cached)
└── reports/                    # Output files (Results)
    ├── webui/                  # UI specific uploads/states
    ├── distances.csv           # Calculated distances between pairs
    └── enrichment_summary.md   # Final statistical report
```

## ⚙️ Configuration & Notes

*   **AlphaFold Structures**: Stored in `alphafold_structures/`. The downloader attempts multiple version suffixes (`v4`, `v3`, etc.) to ensure the best model is found.
*   **Performance**: Distance computation is optimized but can be heavy for datasets with thousands of proteins. 
*   **Custom PTMs**: When uploading custom PTMs, providing "Allowed Residues" (e.g., `K` for Ubiquitination) helps filter out invalid sites automatically.

## 📄 License

[MIT License](LICENSE)
