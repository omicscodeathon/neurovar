# NeuroVar Shiny Application User Guide

## Implementation and  Operation

    - platform-independent
    - all dependencies are installed automatically with the tool
    - Internet requirement: Only if using an online server
    - Third-party tools requirement: R studio or an online Server
    - Pre-requirements: R  installation
    
### Installation

**Step 1 :** Download the application
   
    Download the application folder "NeurVar_Shiny"     
    from [here](https://github.com/omicscodeathon/neurovar/tree/main/NeuroVar_Shiny_Application/NeurVar_Shiny)                                                                                                                      
**Step 2 :** Install dependencies
    
    Install the latest version of R if your computer does not have it: [R](https://cran.r-project.org/bin/windows/base/)   
    
            ### How To Install R and RStudio
            - Go to This Website : https://posit.co/download/rstudio-desktop/ 
![image](https://github.com/omicscodeathon/Exvar/assets/73958439/62b7eda6-c7af-47a2-aec9-fe14aae68e50)

            -  Download R for you specific Operating System (OS). R 4.2.1 version is recommended [Link to version 4.2.1](https://cran.r-project.org/bin/windows/base/old/4.2.1/)

![image](https://github.com/omicscodeathon/Exvar/assets/73958439/258e6366-4cf9-45a9-ba39-ebaf4212af71)

           - Download RStudio

**Step 3 :** Settings

    - Open the application folder (NeurVar_Shiny) using RStudio   

    - Set the folder as a working directory 

    - Unrar the folder "annotation" in the "internal data" folder.

    - Then open the file named "Run_App.R"

**Step 4 :** Run The Application

Run the app code by cicking on the "Run APP" button
 
![image](https://github.com/omicscodeathon/neurovar/assets/73958439/249ad53a-db4c-4bd2-b0d9-9babd98bd8df)

All dependencies will be installed automaticly and the application's dashboard will appear.

**Step 5:**

The tool integrates 3 tabs, each for a different functionality.

In the First tab, select the patient’s disease, next, a list of specific subtypes of the disease with neurological manifestation is provided. More information about each gene is provided including mode of inheritance, description, type, and transcripts. A Link to the official online report validating Gene's association with the disease is also provided.


![Figure 1](https://github.com/omicscodeathon/neurovar/assets/73958439/633a6d98-959f-4522-9f39-d5a3774382e3)


**Step 6:** Visualize the biomarkers expression profile:

   
In the second tab, import the gene expression data file (as a CSV file) and identify the key columns, the log2FC value and p-value to define the differential expression profile. The genes’ expression profiles are summarized in a table and represented in a volcano plot.


![fig2 2](https://github.com/omicscodeathon/neurovar/assets/73958439/917b3134-2d7b-4c37-b529-308d4d5b51c5)

![image](https://github.com/user-attachments/assets/81b82292-6e4f-4eea-8a0b-8c7734dc466b)

**Step 7:** Visualize SNPs and Indels data:

In the third tab,  define the path to the directory containing the genetic variants data files ( as VCF files). Then define the variants type (SNPs or Indels). The VCF files are processed and annotated, and then the variants in the biomarkers are filtered and resumed in a table comparing the reference genome, the control group and the patient group.


![Figure 3 indel](https://github.com/omicscodeathon/neurovar/assets/73958439/a886405a-8282-4624-8f4c-edeedd7f517c)


#### Demonstration Data

Demonstration data are provided [here](https://github.com/omicscodeathon/neurovar/tree/main/demonstration_data)

This data will guide you to organize your data correctly and could be used to test the tool.

## License

Artistic license 2.0
