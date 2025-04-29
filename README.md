# 🧬 SARS-CoV-2 Whole Genome Variant Analysis Report

**Data Source:** Variant annotations from whole genome sequencing of SARS-CoV-2 samples\
**Tools Used:** Bash, Nextflow, R.

## 📘 Overview

This report provides a comprehensive exploration of variants identified from whole genome sequencing (WGS) of SARS-CoV-2 samples. The initial processing and analysis were conducted using the ViralRecon pipeline on Nextflow via Bash scripting, which enabled high-throughput quality control, read mapping, variant calling, and consensus sequence generation. Using visual analytics in R, I aimed to uncover the biological patterns behind mutation distribution, variant types, amino acid substitutions, and more.

------------------------------------------------------------------------

## 1. 🔬 Distribution of Variants per Gene
![variants per gene](https://github.com/user-attachments/assets/319f9d4b-d2ba-4933-ab0d-2f79fb30065c)

### 📌 Interpretation:

-   The **Spike (S) gene** has the highest number of variants. This is significant because the Spike protein mediates entry into host cells and is a primary target of vaccines.
-   Variants are also prominent in **ORF1ab**, which encodes the replicase complex — essential for viral replication.
-   Mutation hotspots here could indicate **adaptive evolution** in response to immune pressure or antiviral therapy.

------------------------------------------------------------------------

## 2. 🧠 Functional Effect of Variants
![proportion of mutation effects](https://github.com/user-attachments/assets/378a27d9-664e-409f-a327-06a838d45a34)

### 📌 Interpretation:

-   **Missense mutations** dominate the dataset, which can lead to amino acid substitutions and potentially alter protein function.
-   A sizable portion are **synonymous**, meaning they do not change the amino acid but may affect translation efficiency.
-   **Stop gained** and **frameshift mutations**, though fewer, could have drastic effects by truncating viral proteins.

------------------------------------------------------------------------

## 3. 🌡️ Heatmap of Mutation Density Across Genome

![heatmap of variant positions by sample](https://github.com/user-attachments/assets/dc8993ad-6af6-42f9-8ad6-b11f32f258b3)

### 📌 Interpretation:

-   Certain genome regions, such as 21,000–25,000 (S gene region), show **higher mutation density**.
-   This suggests evolutionary hotspots and potential selective pressure.
-   Mutational clustering may align with **emerging variants of concern (VOCs)**.

------------------------------------------------------------------------

## 4. 🧬 Top Amino Acid Changes
![amino acid changes](https://github.com/user-attachments/assets/a2776b69-fc02-43b3-a454-820f870220b6)

### 📌 Interpretation:

-   High-frequency mutations such as **D614G**, **E484K**, and **Q677H** were observed.
-   **D614G** increases viral infectivity and has become globally dominant.
-   **E484K** is associated with **immune escape**, often found in Beta and Gamma variants.
-   These mutations warrant close monitoring due to their **functional and epidemiological impact**.

------------------------------------------------------------------------

## 5. 📈 Allele Frequency Distribution
![allele frequency distribution](https://github.com/user-attachments/assets/3f18a5b9-6995-408e-86f5-31f90f3e4201)

### 📌 Interpretation:

-   Most mutations have an **allele frequency close to 1**, indicating they are **fixed or dominant** in those samples.
-   Variants with lower allele frequency may be **emerging** or represent **intra-host diversity**.
-   This insight is crucial for understanding **viral population dynamics** within hosts.

------------------------------------------------------------------------

## 6. 👥 Variants per Sample
![variant count per sample](https://github.com/user-attachments/assets/16340197-af81-4233-9239-871c7bae3b55)

### 📌 Interpretation:

-   Some samples contain more variants, which may be due to **prolonged infection**, **host immune status**, or sequencing depth.
-   Tracking mutation load per sample can help identify **hypermutated viruses** or **host-driven mutation pressure**.

------------------------------------------------------------------------

## 🧩 Summary

This analysis shows that the **Spike protein** is a major site of evolutionary change in SARS-CoV-2. - Mutations like **D614G** and **E484K** have major implications for **transmission and immune escape**. - Variant allele frequencies suggest both **dominant and emerging variants** within the population. This kind of genomic surveillance can guide **public health decisions** and **vaccine updates**.

------------------------------------------------------------------------

## ✍️ Citation

If using this report, please cite as:\
\> [Chijioke Uhegwu]. "SARS-CoV-2 Whole Genome Variant Analysis." GitHub. [https://github.com/ChijiokeUhegwu/viralrecon], 2025.

------------------------------------------------------------------------

## 📬 Contact

For questions or collaboration:\
📧 [chijiokeuhegwu@gmail.com]\
🔗 GitHub: [ChijiokeUhegwu]
