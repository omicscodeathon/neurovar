# NeuroVar Shiny Application User Guide

## Implementation and  Operation

    - platform-independent
    - all dependencies are installed automatically with the tool
    - Internet requirement: Only if using an online server
    - Third-party tools requirement: R studio or an online Server
    - Pre-requirements: R  installation
    
### Installation

**Step 1 :** Download the application
   
    Download the application folder from [here](https://github.com/omicscodeathon/neurovar/tree/main/NeuroVar_Shiny_Application)                                                                                                                      
**Step 2 :** Install dependencies
    
    Install the latest version of R if your computer does not have it: [R](https://cran.r-project.org/bin/windows/base/)   
    
            ### How To Install R and RStudio
            - Go to This Website : https://posit.co/download/rstudio-desktop/ 
![image](https://github.com/omicscodeathon/Exvar/assets/73958439/62b7eda6-c7af-47a2-aec9-fe14aae68e50)

            -  Download R for you specific Operating System (OS). R 4.2.1 version is recommended [Link to version 4.2.1](https://cran.r-project.org/bin/windows/base/old/4.2.1/)

![image](https://github.com/omicscodeathon/Exvar/assets/73958439/258e6366-4cf9-45a9-ba39-ebaf4212af71)

           - Download RStudio

**Step 3 :** Open the application 

    - Open the application folder using RStudio   

    - set the shiny app directory as a working directory   

    - Then open the file named "neurovar.R"

**Step 4 :** Run the app by cicking on the "Run APP" button 
 
![image](https://github.com/omicscodeathon/neurovar/assets/73958439/249ad53a-db4c-4bd2-b0d9-9babd98bd8df)

**Step 5:**

The tool integrates 3 tabs, each for a different functionality.

In the First tab, select the patient’s disease, next, a list of specific subtypes of the disease with neurological manifestation is provided. More information about each gene is provided including mode of inheritance, description, type, and transcripts. A Link to the official online report validating Gene's association with the disease is also provided.

<p align="center">
<img src="https://user-images.githubusercontent.com/73958439/232723944-8e5e658e-bbe5-40e7-92d7-f855ae0400aa.png" width="600" >
</p>

**Step 6:** Visualize the biomarkers expression profile:

   
In teh second tab, import the gene expression data file (as a CSV file) and identify the key columns, the log2FC value and p-value to define the differential expression profile. The genes’ expression profiles are summarized in a table and represented in a volcano plot.

<p align="center">
<img src="https://user-images.githubusercontent.com/73958439/232724064-e2803d44-4381-408e-b0b7-0a9553b8a16b.png" width="600" >
</p>

**Step 7:** Visualize SNPs and Indels data:

In the third tab,  define the path to the directory containing the genetic variants data files ( as VCF files). Then define the variants type (SNPs or Indels). The VCF files are processed and annotated, and then the variants in the biomarkers are filtered and resumed in a table comparing the reference genome, the control group and the patient group.

<p align="center">
<img src="https://user-images.githubusercontent.com/73958439/232724156-3bd91417-89ec-4d1e-a56e-953836a0256b.png" width="600" >
</p>


#### Demonstration Data

Demonstration data are provided [here](https://github.com/omicscodeathon/neurovar/blob/main/demonstration_data.rar)

This data will guide you to organize your data correctly and could be used to test the tool.

## License

Artistic license 2.0
