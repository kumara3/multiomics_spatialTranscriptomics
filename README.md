# Glioma Niches in Spatial Transcriptomics (DMG H3K27M vs GBM)



This repository documents an analysis workflow and key takeaways from a review of publication: Spatial transcriptomics reveals niche-specific enrichment and vulnerabilities of radial glial stem-like cells in malignant gliomas (https://www.nature.com/articles/s41467-023-36707-6). The study is focussed on **diffuse midline glioma (DMG, H3K27M-mutant)** and **glioblastoma (GBM)**, focused on **spatial organization of cellular states**, **niche-specific regulatory programs**, **copy-number variation**, **deconvolution of cell types**, and **isoform/splicing diversity** across glioma niches. :contentReference[oaicite:0]{index=0}

---

## Project goals

- Integrate **spatial transcriptomics** from multiple glioma samples in both **transcriptional** and **spatial** spaces. :contentReference[oaicite:1]{index=1}  
- Identify **clinically relevant transcriptional programs**, **multi-cellular ecosystems**, and **niche-defined modules** across tumor regions. :contentReference[oaicite:2]{index=2}  
- Use **CNV inference** and **cell-type deconvolution** to distinguish malignant vs microenvironmental content and map it spatially. :contentReference[oaicite:3]{index=3}  

---

## Background (from the PDF)

- **DMG**: typically pediatric, midline brain regions; associated with **H3K27M**, often alongside **TP53 mutation** and **PDGFRA amplification**. :contentReference[oaicite:4]{index=4}  
- **GBM**: typically adult, often frontal/temporal; driven by complex genetic events. :contentReference[oaicite:5]{index=5}  
- Both have poor prognosis under standard treatment paradigms (noted in the deck). :contentReference[oaicite:6]{index=6}  

---

## Data summary

- **11 total samples** were integrated (DMG + GBM; includes a peritumor tissue sample for GBM5_2). :contentReference[oaicite:7]{index=7}  
- The workflow references both **2nd-gen short-read** and **3rd-gen long-read sequencing**. :contentReference[oaicite:8]{index=8}  

---

## Methods overview

### 1) Preprocessing and clustering
- Input: Cell Ranger filtered raw counts
- Primary workflow: **Seurat** across 11 samples
- Integration:
  - **Harmony** for transcriptional space integration
  - **BANKSY** for spatial integration :contentReference[oaicite:9]{index=9}  

### 2) Gene set scoring → transcriptional programs
- Differential expression marker sets were built per cluster/cell type.
- A subset of gene sets (e.g., >30 genes) were scored per cell via **Seurat `AddModuleScore`**.
- Modules were clustered (Pearson correlation distance; hierarchical clustering) and “forced” into 4 clusters via `cutree`. :contentReference[oaicite:10]{index=10}  

### 3) CNV inference (malignant vs reference)
- **inferCNV** used after integration.
- Reference set described as cells expressing markers such as **ZIC1, MBP, GABRA1** in specific clusters; all other cells treated as tumor. :contentReference[oaicite:11]{index=11}  

### 4) Niche/module interpretation
- The deck maps module clusters to niche-like labels such as:
  - **Tumor core**
  - **Vascular**
  - **Invasive**
  - (Summary also mentions hypoxic among niches) :contentReference[oaicite:12]{index=12}  

### 5) Cell-type deconvolution (spatial)
- Tool: **RCTD** (Robust decomposition of cell type mixtures in spatial transcriptomics).
- Reference datasets cited in the deck include:
  - Bhaduri et al. (GBM scRNA-seq)
  - Filbin et al. (DIPG / H3K27M scRNA-seq)
  - Aldinger et al. (human cerebellum snRNA-seq)
  - Nowakowski et al. (human cortex scRNA-seq) :contentReference[oaicite:13]{index=13}  

---

## Notable findings

- Transcriptional programs include cell-cycle, astrocytic differentiation (AC-like), oligodendrocytic differentiation (OC-like), and OPC-like programs (e.g., PDGFRA/CSPG4). :contentReference[oaicite:14]{index=14}  
- The study describes **RG-like (radial glial-like)** signatures being identified/defined using RG markers within astrocyte-like cells (AC2) and notes marker examples: **TNC, HOPX, PTPRZ1, VIM**. :contentReference[oaicite:15]{index=15}  
- Deconvolution suggests vascular and immune components (e.g., endothelial/pericytes, microglia) appear across transcriptional programs. :contentReference[oaicite:16]{index=16}  

---

## Known issues / caveats
- “SpotClean and BayesTME error” is explicitly noted. :contentReference[oaicite:17]{index=17}  
- A stated drawback: inclusion of samples spanning **peritumor and normal brain areas** with comparable anatomic structures may complicate interpretation. :contentReference[oaicite:18]{index=18}  

---


