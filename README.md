# Glioma Niches in Spatial Transcriptomics (DMG H3K27M vs GBM)

This repository documents an analysis workflow and key takeaways from a review of publication: Spatial transcriptomics reveals niche-specific enrichment and vulnerabilities of radial glial stem-like cells in malignant gliomas (https://www.nature.com/articles/s41467-023-36707-6). The study is focussed on **diffuse midline glioma (DMG, H3K27M-mutant)** and **glioblastoma (GBM)**, focused on **spatial organization of cellular states**, **niche-specific regulatory programs**, **copy-number variation**, **deconvolution of cell types**, and **isoform/splicing diversity** across glioma niches. :contentReference[oaicite:0]{index=0}

---

## Project goals

- Integrate **spatial transcriptomics** from multiple glioma samples in both **transcriptional** and **spatial** spaces.   
- Identify **clinically relevant transcriptional programs**, **multi-cellular ecosystems**, and **niche-defined modules** across tumor regions. 
- Use **CNV inference** and **cell-type deconvolution** to distinguish malignant vs microenvironmental content and map it spatially.  

---

## Background 

- **DMG**: typically pediatric, midline brain regions; associated with **H3K27M**, often alongside **TP53 mutation** and **PDGFRA amplification**. 
- **GBM**: typically adult, often frontal/temporal; driven by complex genetic events.  
- Both have poor prognosis under standard treatment paradigms (noted in the deck).

---

## Data summary

- **11 total samples** were integrated (DMG + GBM; includes a peritumor tissue sample for GBM5_2). 
- The workflow references both **2nd-gen short-read** and **3rd-gen long-read sequencing**.

---

## Methods overview

### 1) Preprocessing and clustering
- Input: Cell Ranger filtered raw counts
- Primary workflow: **Seurat** across 11 samples
- Integration:
  - **Harmony** for transcriptional space integration
  - **BANKSY** for spatial integration 

### 2) Gene set scoring → transcriptional programs
- Differential expression marker sets were built per cluster/cell type.
- A subset of gene sets (e.g., >30 genes) were scored per cell via **Seurat `AddModuleScore`**.
- Modules were clustered (Pearson correlation distance; hierarchical clustering) and “forced” into 4 clusters via `cutree`.

### 3) CNV inference (malignant vs reference)
- **inferCNV** used after integration.
- Reference set described as cells expressing markers such as **ZIC1, MBP, GABRA1** in specific clusters; all other cells treated as tumor. 

### 4) Niche/module interpretation
- The deck maps module clusters to niche-like labels such as:
  - **Tumor core**
  - **Vascular**
  - **Invasive**

### 5) Cell-type deconvolution (spatial)
- Tool: **RCTD** (Robust decomposition of cell type mixtures in spatial transcriptomics).
- Reference datasets cited in the deck include:
  - Bhaduri et al. (GBM scRNA-seq)
  - Filbin et al. (DIPG / H3K27M scRNA-seq)
  - Aldinger et al. (human cerebellum snRNA-seq)
  - Nowakowski et al. (human cortex scRNA-seq)
---

## Notable findings

- Transcriptional programs include cell-cycle, astrocytic differentiation (AC-like), oligodendrocytic differentiation (OC-like), and OPC-like programs (e.g., PDGFRA/CSPG4).   
- The study describes **RG-like (radial glial-like)** signatures being identified/defined using RG markers within astrocyte-like cells (AC2) and notes marker examples: **TNC, HOPX, PTPRZ1, VIM**.  
- Deconvolution suggests vascular and immune components (e.g., endothelial/pericytes, microglia) appear across transcriptional programs.  

---

## Known issues / caveats
- “SpotClean and BayesTME error” is explicitly noted. 
- A stated drawback: inclusion of samples spanning **peritumor and normal brain areas** with comparable anatomic structures may complicate interpretation.
---


