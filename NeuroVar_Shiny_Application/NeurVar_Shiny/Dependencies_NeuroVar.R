#requirement
if (!require("shiny")) install.packages("shiny")
if (!require("dplyr")) install.packages("dplyr")
if (!require("readr")) install.packages("readr")
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("purrr")) install.packages("purrr")
if (!require("vcfR")) install.packages("vcfR")
if (!require("bslib")) install.packages("bslib")
if (!require("stringr")) install.packages("stringr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("shinydashboard")) install.packages("shinydashboard")
if (!require("DT")) install.packages("DT")
if (!require("vcfR")) install.packages("vcfR")
if (!require("sqldf")) install.packages("sqldf")
if (!require("fs")) install.packages("fs")
library(shiny)
library(dplyr)
library(readr)
library(tidyverse)
library(purrr)
library(vcfR)
library(bslib)
library(stringr)
library(ggplot2)
library(shinydashboard) 
library(DT)
library(sqldf)
library(data.table)
library(fs)
full_list <- read_csv("internal_data/full_list.csv")
unzip_dir <- "internal_data"
zip_file1 <- "internal_data/annotation_part1.zip"
unzip(zip_file1, exdir = unzip_dir)
zip_file2 <- "internal_data/annotation_part2.zip"
unzip(zip_file2, exdir = unzip_dir)
annotation1_read <- read.table("internal_data/annotation_part1.txt", sep = "\t", header = TRUE)
annotation2_read <- read.table("internal_data/annotation_part2.txt", sep = "\t", header = TRUE)
annotation<- rbind(annotation1_read, annotation2_read)
