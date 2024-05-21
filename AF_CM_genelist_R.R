getwd()
setwd("C:/Users/shsan/Documents/Medical School/TTC and Clerkships (2024-2025)/Research Elective Summer 2024")

library(dplyr)
library(readxl)
library(writexl)

#base text file pulled from NCBI using (atrial fibrillation[All Fields] AND cardiomyopathy[All Fields]) AND ("genetype protein coding"[Properties] AND alive[prop])
afcmgenelist1 <- read.delim("~/Medical School/TTC and Clerkships (2024-2025)/Research Elective Summer 2024/AF_CM_genelist_NCBI.txt")
afcmgenelist1<- data.frame(afcmgenelist1)
View(afcmgenelist1)
afcmgenelist1$descriptors<- paste(afcmgenelist1$description, afcmgenelist1$other_designations, sep = " ")

afcmgenelist1$Classification <- 
  ifelse(grepl("channel|voltage", afcmgenelist1$descriptors),afcmgenelist1$Classification <- 1,
       ifelse (grepl("muscle|desmin|lamin|myosin|sarcomer|titin", afcmgenelist1$descriptors), afcmgenelist1$Classification <- 2,0)
)

list1<- afcmgenelist1[c('Symbol','Classification')]

#Used dbGap -> atrial fibrillation -> phenotype datasets. 
#One study was the Vanderbilt study, so pulled panel of genes from the Yoneda et al. 2021 paper.
afcmgenelist2 <- read.table("~/Medical School/TTC and Clerkships (2024-2025)/Research Elective Summer 2024/AF_CM_genelist.txt", quote="\"", comment.char="")
list2<- unlist(afcmgenelist2)
list2<- data.frame(list2)
View(list2)

write_xlsx(list1, path = "AF_CM_Genelist1_final.xlsx", col_names = FALSE, format_headers = FALSE)
write_xlsx(list2, path = "AF_CM_Genelist2_final.xlsx", col_names = FALSE, format_headers = FALSE)

#merging excel files into a combined list after manual review and updated excel sheets.
#Cells may contain notes/annotations.
list1_final <- read_excel("~/Medical School/TTC and Clerkships (2024-2025)/Research Elective Summer 2024/AF_CM_Genelist1_final.xlsx")
list2_final <- read_excel("~/Medical School/TTC and Clerkships (2024-2025)/Research Elective Summer 2024/AF_CM_Genelist2_final.xlsx")

comb_list <- rbind(list1_final, list2_final[-1,])
comb_list<- unique(comb_list)
af_cm_genelist_combined<- write_xlsx(comb_list, "C:/Users/shsan/Documents/Medical School/TTC and Clerkships (2024-2025)/Research Elective Summer 2024/AF_CM_Genelist_1+2_combined.xlsx")
sarcomeric_count<- sum(comb_list$Classification == 2)
ion_channel_count <- sum(comb_list$Classification == 1)
non_ion_channel_sarcomeric_count <- sum(comb_list$Classification == 0)

