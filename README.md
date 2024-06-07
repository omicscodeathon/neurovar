# NeuroVar: A Genetic Expression and Variation data visualization tool for Neurological diseases’ biomarkers

## Table of Contents
1. [Background](#Background)
2. [About NeuroVar](#About-NeuroVar)
3. [Implementation and  Operation](#Implementation-and-Operation)
4. [Installation](#Installation)
5. [Usage](#Usage)
6. [Demonstration](#Demonstration)
7. [Citation](#Citation)
8. [License](#License)
9. [Contributors](#Contributors)
<br>

## Background

The expanding availability of large-scale genomic data and the growing interest in uncovering gene-disease associations call for efficient tools to visualize and evaluate gene expression and genetic variation data.  

## About NeuroVar
NeuroVar is a novel tool for visualizing genetic variation (Single nucleotide polymorphisms and insertions/deletions) and gene expression data related to neurological diseases. The tool is available as a desktop application that does not require any computational skills to use.

The current version includes :

* 11 Neurological Syndromes are integrated, including:

              ★ Epilepsy 
              ★ Amyotrophic lateral sclerosis 
              ★ Intellectual disability and Autism spectrum disorder 
              ★ Brain malformation syndrome 
              ★ Syndromic Disorders 
              ★ Cerebral Palsy 
              ★ RASopathy 
              ★ Aminoacidopathy 
              ★ CranioFacial Malformations 
              ★ Parkinson disease 
              ★ PHARC syndrome 

* and 7 Non-Neurological diseases with neurological manifestations are integrated, including:

              ★ Peroxisomal disorders  
              ★ Hereditary cancer 
              ★ Mitochondrial disease 
              ★ Retina-related disorders 
              ★ General Gene Curation
              ★ Hearing Loss 
              ★ Fatty Acid Oxidation Disorders


## Implementation and  Operation

NeuroVar is available as a Shiny Application and a desktop application.

1. Shiny application

Guide :  
 
2. Desktop application

Guide :  [Click Here](https://github.com/omicscodeathon/neurovar/tree/main/NeuroVar_Desktop_Application)      
    
## Inputs 

NeuroVar requires gene expression CSV files and genetic variants VCF files as inputs.


## Usage

The tool integrates 3 tabs, each for a different functionality,  as explained below:

1. First the user selects the patient’s disease, next, a list of specific subtypes of the disease with neurological manifestation is provided. More information about each gene is provided including mode of inheritance, description, type, and transcripts. A Link to the official online report validating Gene's association with the disease is also provided.

<p align="center">
<img src="https://user-images.githubusercontent.com/73958439/232723944-8e5e658e-bbe5-40e7-92d7-f855ae0400aa.png" width="600" >
</p>

2. Visualize the biomarkers expression profile:
After importing a CSV file and identifying the key columns, the log2FC value and p-value to define the differential expression profile are requested. The genes’ expression profiles are summarized in a table and represented in a volcano plot.

<p align="center">
<img src="https://user-images.githubusercontent.com/73958439/232724064-e2803d44-4381-408e-b0b7-0a9553b8a16b.png" width="600" >
</p>

3. Visualize SNPs and Indels data:
The user is requested to define the path to the directory containing the VCF files. The user needs to define the variants type (SNPs or Indels). The VCF files are processed and annotated, and then the variants in the biomarkers are filtered and resumed in a table comparing the reference genome, the control group and the patient group.

<p align="center">
<img src="https://user-images.githubusercontent.com/73958439/232724156-3bd91417-89ec-4d1e-a56e-953836a0256b.png" width="600" >
</p>


## Demonstration Data

Demonstration data are provided [here](https://github.com/omicscodeathon/neurovar/blob/main/demonstration_data.rar)

This data will guide you to organize your data correctly and could be used to test the tool.

## Citation

Hiba Ben Aribi, Najla Abassi and Olaitan I. Awe. NeuroVar: A Genetic Expression and Variation data visualization tool for Neurological diseases’ biomarkers https://github.com/omicscodeathon/neurovar



## Contributors

Hiba Ben Aribi,  Faculty of Sciences of Tunis, University of Tunis El Manar, Tunis, Tunisia  (Developer - Writer).

Najla Abassi, Pasteur Institute of Tunis,, Tunisia  (Developer - Writer).

Olaitan I. Awe, African Society for Bioinformatics and Computational Biology (ASBCB), South Africa.
