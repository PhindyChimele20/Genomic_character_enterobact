

#  README

## 1. Data sources

This task is based on on publicly available sequencing data from the study "Genomic characterization of Enterobacter sp. PGRG2 and Achromobacter insolitus PGRG5: bacterial strains isolated from soil present near electronics manufacture industry for heavy metal remediation" which focused on understanding the gene expression of Saccharolobus islandicus relative to its different cell-cycle stages. Cultures of S. islandicus were made to begin the cell cycle at once by a 6-hour treatment with acetic acid. Whole genome sequencing was done on both extracted strains where a sequencing library with paired-end reads (2 × 150 bp) was created and subsequently subjected to sequencing on the NovaSeq 6000 platform provided by Illumina (3).
The subsampled FASTQs are used as the inputs for the workflow.

---

## 2. How to download

The samples were obtained by first fetching the sample’s SRA record (SRX) on GEO using geofetch, then downloaded the fastq of the selected samples using their corresponding SRR accession.
### Code for downloading

```bash
geofetch -i GSE296035 --just-metadata
SRRS=("SRX28623476" "	SRX28623471" "SRX28623484")

for SRR in "${SRRS[@]}"; do
    echo "Downloading $SRR ..."
    prefetch "$SRR"
    fastq-dump --gzip --split-files "$SRR"
done
```


---

## 3. Pre-processing 

From the GEO where the raw data fastq files are, we selected 3 samples and downloaded using their SRR accessions as per above script (Code for Downloaded)

1. **STEP 1** ...

Example:

```bash
CODE TO SUBSAMPLE
```


---

## 4. How the workflow works
The workflow files is stored in workflow/ and it is divided into different steps:
DESCRIBE THE WORKFLOW HERE - NOTE THE BELOW ARE JUST EXAMPLES, REPLACE WITH YOUR OWN - YOURS CAN TAKE A VERY DIFFERENT FORMAT
The workflow files is stored in `workflow/`.

---

### Step 1 – Alignment stage

**Purpose:** The workflow takes each FASTQ file, maps reads to the reference using Bowtie2, converts the alignments to BAM format with Samtools, and deletes the intermediate SAM files.
**Tools:** `Bowtie2`, `Samtools`
**Inputs:** Subsampled FASTQ files (from `data/fastq_subsampled/`)
**Outputs:** bam and bai files, QC reports (`.html`, `.txt`)
**Command:**

```bash
for i in "${samples[@]}"; do
    bowtie2 -x ref_index -U "${i}.fastq" -S "${i}.sam"
    samtools view -b "${i}.sam" > "${i}.bam"
    rm "${i}.sam"
done

```

---

### Step 2 - Quantification

**Purpose:** This part of the workflow takes all the aligned reads (BAM files), uses the reference annotation (sequence.gtf), and produces a matrix of raw read counts per gene per sample.The output file is used fow downstream differential expression analysis
**Tools:** 'featureCounts'
**Inputs:** bam files
**Outputs:** counts matrix (.txt)
**Command:**

```bash
featureCounts -a sequence.gtf -F GTF -o counts_2.txt -T 14 *.bam

```
---

### Step 3 – Differential gene expression analysis

**Purpose:** the pipeline loads is done using in R where the raw gene counts are loaded, filter low expression genes, normalizes counts (TMM), calculates rough log2 fold changes, and uses limma-voom to perform proper differential expression testing with statistics (p-values, adjusted p-values)
**Tools:** 'edgeR','limma'
**Inputs:** count matrix
**Outputs:** Normalized expression values (CPM), Log2 fold-changes (manual + statistical), differential expression test results ('fit2')
**Command:**
```bash
#installing edgeR library
BiocManager::install("edgeR")
library(edgeR)#importing library

#importing the data, making Geneid the first col, making header true
counts = read.table("counts.txt",row.names = 1, header=T)
counts<-counts[, -(1:5)]# removing irrelevant columns
dge <- DGEList(counts = counts)

# Convert to CPM
cpm_counts <- cpm(dge)

# Keep genes with at least 1 CPM for all samples
keep <- rowSums(cpm_counts >= 1) == 3
dge <- dge[keep, , keep.lib.sizes=FALSE]

# TMM normalisation
dge <- calcNormFactors(dge, method="TMM")
###############################################
# Define columns for each group (updating sample names from metadata)
# Rename columns of counts
colnames(counts) <- c("M_G1", "S", "G2")
colnames(counts)

norm_counts <- cpm(dge, normalized.lib.sizes = TRUE)

# Mean log2 expression per group
mean_MG1 <- rowMeans(log2(norm_counts[, M_G1,drop = FALSE] + 1))
mean_S   <- rowMeans(log2(norm_counts[, S,drop = FALSE] + 1))
mean_G2  <- rowMeans(log2(norm_counts[, G2,drop = FALSE] + 1))

# Calculate log2 fold changes
log2FC_MG1_vs_S  <- mean_MG1 - mean_S
log2FC_S_vs_G2   <- mean_S   - mean_G2
log2FC_G2_vs_MG1 <- mean_G2  - mean_MG1

# Combine into one results table, keeping the same index as norm_counts
results <- data.frame(
  log2FC_MG1_vs_S  = log2FC_MG1_vs_S,
  log2FC_S_vs_G2   = log2FC_S_vs_G2,
  log2FC_G2_vs_MG1 = log2FC_G2_vs_MG1
)

# Preserve the gene IDs as rownames
rownames(results) <- rownames(norm_counts)

head(results)
##############
library(limma)

# --- Define sample groups ---
# Replace these with the actual sample labels in your counts table
group <- factor(c(rep("M_G1", 1), rep( "S", 1), rep("G2", 1)))

# --- Design matrix ---
design<- model.matrix(~0 + group)
colnames(design) <- levels(group)

# --- Transform counts with voom ---
v <- voom(dge, design, plot = FALSE)

# --- Fit linear model ---
fit <- lmFit(v, design)

# --- Define contrasts ---
contrasts <- makeContrasts(
  M_G1_vs_S  = M_G1 - S,
  S_vs_G2    = S - G2,
  G2_vs_M_G1 = G2 - M_G1,
  levels = design
)

fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2) # No residual degrees of freedom in linear model fits


```
