# NCATS ASMS 384-Well Pooling Engine

An interactive web application designed to automate the design of pooled compound libraries for Affinity Selection-Mass Spectrometry (AS-MS) screening assays. This tool ingests structural library data, predicts mass spectrometry ionization behavior, and generates optimized 384-well plate layouts that maximize mass diversity while preventing isobaric clashing.

## Key Features

* **Dynamic Ionization Prediction (SMARTS):** Automatically predicts whether a compound will ionize in positive or negative mode using RDKit-driven SMILES Arbitrary Target Specifications (SMARTS) pattern matching. 
    * *Positive Mode Rule:* Evaluates basic, non-amide, non-aromatic nitrogens utilizing the SMARTS query: `[NX3;H2,H1,H0;!$(NC=O);!$(N-[#6a])]` (protonates to [M+H]+ with target m/z = Exact Mass + 1.0073).
    * *Negative Mode Rule:* Evaluates acidic carboxylic or sulfonic groups utilizing the SMARTS query: `[C,S](=[O,S])[O;H1,-1]` (deprotonates to [M-H]- with target m/z = Exact Mass - 1.0073).
    * *Amphoteric & Neutral Fallback:* Automatically defaults zwitterionic, amphoteric, or neutral compounds to positive mode to align with standard screening protocols.
* **Mass-Diversity Stride Algorithm:** Resolves the critical bottleneck of isobaric overlapping. The engine sorts compounds by ascending m/z within each polarity group and applies a modular stride distribution to space out similar masses across different wells, maximizing the delta-m/z within every single pool.
* **Continuous 384-Well Plate Mapping:** Consolidates both ionization pools back-to-back onto a single physical plate footprint (A01 to P24). Placing all negative pools sequentially first followed immediately by positive pools allows downstream mass spectrometers to run dedicated, uninterrupted polarity methods blocks, preventing the instrument from having to continually switch polarities.
* **Robust Chemistry Parsing & Sanitization:** Accepts Spotfire .SDF library exports. Automatically extracts structural properties (checking fields like `SAMPLE_ID`, `ID`, `Name`, and `SMILES`), systematically strips trailing lot/batch suffixes (e.g., `-01`, `-07`) at the root loop to enable clean downstream database integration, and filters out corrupted "dead" records, salt fragments, or molecules missing structural data.
* **Self-Contained Interactive Visual Layouts:** Exports a standalone HTML visual plate map with an embedded JSON payload that requires zero server architecture. Color-coded by predicted polarity (pink for negative, blue for positive), the map allows scientists to click any well to inspect exact masses, target m/z values, and clean vector-based SVG chemical structures rendered via a segfault-safe C++ drawing engine.
* **Automation-Ready Manifests:** Generates structured, GLP-compliant layout CSV manifests matching exact plate coordination paths, designed to serve as direct layouts for CoMa database ingestion or acoustic liquid handlers (such as the Beckman Coulter Echo).

## Pooling & Stride Mathematics

To prevent mass-clashing, the algorithm distributes compounds using modular arithmetic. Let C be the total number of compounds in a given ionization mode, N be the fixed pool size (10 compounds per well), and i represent the sorted index of a compound (0, 1, 2, ..., C-1). 

The total number of required pools/wells (M) is calculated via:

M = floor(C / N)

For any single compound index i, its physical well assignment and internal pool sub-index (layer 1 to 10) are determined by:

Well Index = i mod M

Well Sub-Index = floor(i / M) + 1

This mathematical spacing ensures that the compounds sharing any individual well are separated by a minimum stride of M sorted steps across the weight distribution of the library.

## Technical Stack

* **Language:** Python 3.10+
* **Framework:** Streamlit
* **Cheminformatics Core:** RDKit
* **Data Handling:** Pandas
* **Deployment:** GitHub / Streamlit Community Cloud

## Background

Developed to eliminate manual calculations and complex plate-mapping errors in high-throughput drug screening workflows. Specifically engineered to support AS-MS screening validation of cherry-picked compound libraries against target proteins (such as His-tagged proteins bound to high-affinity Ni-NTA magnetic beads), this engine bridges the gap between raw chemical databases and automated, nanoliter-dispensing liquid handlers. By automating structural parsing and physical mapping, this tool transforms chaotic library exports into pristine, mass-spectrometer-optimized screening plates in a single click.

---

Made by Colin Gordy