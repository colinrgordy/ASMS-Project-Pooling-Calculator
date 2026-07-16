# NCATS ASMS 384-Well Pooling Engine

An interactive web application designed to automate the design of pooled compound libraries for Affinity Selection-Mass Spectrometry (AS-MS) screening assays[span_4](start_span)[span_4](end_span). This tool ingests structural library data, predicts mass spectrometry ionization behavior, and generates optimized 384-well plate layouts that maximize mass diversity while preventing isobaric clashing[span_5](start_span)[span_5](end_span).

## Key Features

* **Dynamic Ionization Prediction (SMARTS):** Automatically predicts whether a compound will ionize in positive or negative mode using RDKit-driven SMILES Arbitrary Target Specifications (SMARTS) pattern matching[span_6](start_span)[span_6](end_span). 
    * *Positive Mode Rule:* Evaluates basic, non-amide, non-aromatic nitrogens utilizing the SMARTS query: `[NX3;H2,H1,H0;!$(NC=O);!$(N-[#6a])]` (protonates to $[M+H]^+$ with target $m/z = \text{Exact Mass} + 1.0073$)[span_7](start_span)[span_7](end_span).
    * *Negative Mode Rule:* Evaluates acidic carboxylic or sulfonic groups utilizing the SMARTS query: `[C,S](=[O,S])[O;H1,-1]` (deprotonates to $[M-H]^-$ with target $m/z = \text{Exact Mass} - 1.0073$)[span_8](start_span)[span_8](end_span).
    * *Amphoteric & Neutral Fallback:* Automatically defaults zwitterionic, amphoteric, or neutral compounds to positive mode to align with standard screening protocols[span_9](start_span)[span_9](end_span).
* **Mass-Diversity Stride Algorithm:** Resolves the critical bottleneck of isobaric overlapping[span_10](start_span)[span_10](end_span). The engine sorts compounds by ascending $m/z$ within each polarity group and applies a modular stride distribution to space out similar masses across different wells, maximizing the delta-$m/z$ within every single pool[span_11](start_span)[span_11](end_span).
* **Continuous 384-Well Plate Mapping:** Consolidates both ionization pools back-to-back onto a single physical plate footprint (A01 to P24)[span_12](start_span)[span_12](end_span). Placing all negative pools sequentially first followed immediately by positive pools allows downstream mass spectrometers to run dedicated, uninterrupted polarity methods blocks, preventing the instrument from having to continually switch polarities[span_13](start_span)[span_13](end_span).
* **Robust Chemistry Parsing & Sanitization:** Accepts Spotfire `.sdf` library exports[span_14](start_span)[span_14](end_span). Automatically extracts structural properties (checking fields like `SAMPLE_ID`, `ID`, `Name`, and `SMILES`), systematically strips trailing lot/batch suffixes (e.g., `-01`, `-07`) at the root loop to enable clean downstream database integration, and filters out corrupted "dead" records, salt fragments, or molecules missing structural data[span_15](start_span)[span_15](end_span).
* **Self-Contained Interactive Visual Layouts:** Exports a standalone HTML visual plate map with an embedded JSON payload that requires zero server architecture[span_16](start_span)[span_16](end_span). Color-coded by predicted polarity (pink for negative, blue for positive), the map allows scientists to click any well to inspect exact masses, target $m/z$ values, and clean vector-based SVG chemical structures rendered via a segfault-safe C++ drawing engine[span_17](start_span)[span_17](end_span).
* **Automation-Ready Manifests:** Generates structured, GLP-compliant layout CSV manifests matching exact plate coordination paths, designed to serve as direct layouts for CoMa database ingestion or acoustic liquid handlers (such as the Beckman Coulter Echo)[span_18](start_span)[span_18](end_span).

## Pooling & Stride Mathematics

To prevent mass-clashing, the algorithm distributes compounds using modular arithmetic[span_19](start_span)[span_19](end_span). Let $C$ be the total number of compounds in a given ionization mode, $N$ be the fixed pool size (10 compounds per well), and $i$ represent the sorted index of a compound ($0, 1, 2, \dots, C-1$)[span_20](start_span)[span_20](end_span). 

The total number of required pools/wells ($M$) is calculated via:

$$M = \left\lfloor \frac{C}{N} \right\rfloor$$

For any single compound index $i$, its physical well assignment and internal pool sub-index (layer 1 to 10) are determined by[span_21](start_span)[span_21](end_span):

$$\text{Well Index} = i \pmod M$$

$$\text{Well Sub-Index} = \left\lfloor \frac{i}{M} \right\rfloor + 1$$

This mathematical spacing ensures that the compounds sharing any individual well are separated by a minimum stride of $M$ sorted steps across the weight distribution of the library[span_22](start_span)[span_22](end_span).

## Technical Stack

* **Language:** Python 3.10+[span_23](start_span)[span_23](end_span)
* **Framework:** Streamlit[span_24](start_span)[span_24](end_span)
* **Cheminformatics Core:** RDKit[span_25](start_span)[span_25](end_span)
* **Data Handling:** Pandas[span_26](start_span)[span_26](end_span)
* **Deployment:** GitHub / Streamlit Community Cloud[span_27](start_span)[span_27](end_span)

## Background

Developed to eliminate manual calculations and complex plate-mapping errors in high-throughput drug screening workflows[span_28](start_span)[span_28](end_span). Specifically engineered to support AS-MS screening validation of cherry-picked compound libraries against target proteins (such as His-tagged proteins bound to high-affinity Ni-NTA magnetic beads), this engine bridges the gap between raw chemical databases and automated, nanoliter-dispensing liquid handlers[span_29](start_span)[span_29](end_span). By automating structural parsing and physical mapping, this tool transforms chaotic library exports into pristine, mass-spectrometer-optimized screening plates in a single click[span_30](start_span)[span_30](end_span).

---

Made by Colin Gordy [span_31](start_span)[span_31](end_span)