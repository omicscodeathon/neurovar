# NeuroVar: A Genetic Expression and Variation data visualization tool for Neurological diseases’ biomarkers

## Background

The expanding availability of large-scale genomic data and the growing interest in uncovering gene-disease associations call for efficient tools to visualize and evaluate gene expression and genetic variation data.  

## About NeuroVar
NeuroVar is a novel tool for visualizing genetic variation (Single nucleotide polymorphisms and insertions/deletions) and gene expression data related to neurological diseases. The tool is available as a desktop application that does not require any computational skills to use.

The current version includes multiple neurological disorders  and diseases with neurological manifestations, including:

★	Epilepsy   ★	amyotrophic lateral sclerosis  ★	autism spectrum disorder   ★	intellectual disability  ★	brain malformation syndrome  ★	Syndromic Disorders   ★	Cerebral Palsy   ★	RASopathy    ★	Amnioacidopathy   ★	Peroxisomal disorders   ★	Hereditary cancer   ★	Mitochondrial disease   ★	Retina related disorders  ★	Hearing LossNeuropathy   ★	neurodevelopmental disorder   ★	fragile X syndrome  ★	frontotemporal dementia   ★	glycine encephalopathy   ★	Griscelli syndrome   ★	Joubert syndrome    ★	Microcephaly   ★	neurodevelopmental disorder   ★	neurofibromatosis   ★	Parkinson disease   ★	PHARC syndrome    ★	tuberous sclerosis  ★	neuroblastoma.## Tool development methodology


## Tool development methodology
The genes associated with the disease are downloaded from the ClinGen database (39). Multiple neurological disorders and diseases with neurological manifestations were selected.
The pipeline is developed using the R language and converted into a shiny application using multiple R packages.
The shiny app was then transformed into a cross-platform desktop application using useing Electron and Node.js.

## Result

1. First the user select the patient’s disease, next, a list of specific subtypes of the disease with neurological manifestation is provided. More information about each gene is provided including mode of inheritance, description, type, and transcripts. A links for the official online report validating the Gene's association with the disease is also provided. 

![image](https://user-images.githubusercontent.com/73958439/232723944-8e5e658e-bbe5-40e7-92d7-f855ae0400aa.png)


3.Visualize the biomarkers expression profile: 
After importing a CSV file and identifying the key columns, the log2FC value and p-value to define the differential expression profile are requested. The genes’ expression profiles are summarized in a table and represented in a volcano plot.

![image](https://user-images.githubusercontent.com/73958439/232724064-e2803d44-4381-408e-b0b7-0a9553b8a16b.png)


3. Visualize SNPs and Indels data:
The user is requested to define the path to the directory containing the VCF files. The user needs to define the variants type (SNPs or Indels). The VCF files are processed and annotated, and then the variants in the biomarkers are filtered and resumed in a table comparing the reference genome, the controls group and the patients group.

![image](https://user-images.githubusercontent.com/73958439/232724156-3bd91417-89ec-4d1e-a56e-953836a0256b.png)



## Implementation and  Operation

The software is platform-independent. All requirements are installed automatically. It requires gene expression CSV files and genetic variants VCF files as inputs.
Note : The files are expected to be divided into two folders, named ''controls'' and "patients", containing the VCF files of the controls and patients groups respectively.

## Pipeline Validation / Case study

To validate the pipeline, a case study was performed on the public dataset SRP149638 (36) available on the SRA database (37). The file’s preprocessing, genetic expression analysis, and variant calling were performed using the Exvar R package (38). 
The dataset corresponds to RNA sequencing data (Expression profiling by high throughput sequencing) from the peripheral blood mononuclear cells (PBMC) from healthy donors and Amyotrophic Lateral Sclerosis (ALS) patients.
The expression analysis shows that only one ALS biomarker is differentially expressed which is TUBA4A gene, a biomarker of ALS type 22. TUBA4A is a protein-coding gene inherited in autosomal dominant mode.

## License: 
Artistic license 2.0 

## Contributors

Hiba Ben Aribi, UTM, Tunisia  (Developer - Writer).

Najla Abassi, UTM, Tunisia  (Developer - Writer).

Olaitan I. Awe, ASBCB, South Africa.
