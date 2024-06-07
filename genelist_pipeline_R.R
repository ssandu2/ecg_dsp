getwd()
setwd("C:/Users/shsan/Documents/Medical School/Research Elective Summer 2024")

library(dplyr)
library(readxl)
library(writexl)
library(httr2)
library(jsonlite)

#AUTOMATING GENE COLLLECTION
#functions--------------------------------------
#obtaining gene lists gen
search_genes <- function (term, retmax = 100){
  url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?",
                "db=gene&term=", URLencode(term), "&retmax=", retmax, "&retmode=json")
  
  search_response<- request(url) %>% req_perform()
  json_response <- resp_body_string(search_response) %>% fromJSON()
  gene_ids <- json_response$esearchresult$idlist
  return(gene_ids)
}

obtain<- function(gene_ids){
  ids<- paste(gene_ids, collapse = ",")
  url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?", "db=gene&id=",
                ids, "&retmode=json")
  
  search_response2<- request(url) %>% req_perform() %>% resp_body_string()
  json_response2 <- fromJSON(search_response2)
  gene_info_prep <- json_response2$result
  gene_info_prep <- gene_info_prep[!names(gene_info_prep) %in% "uids"]
  return(gene_info_prep)
}
#extracting gene info
summarize <- function(summary) {
  list(gene_name <- summary$name,
       gene_aliases <- summary$otheraliases,
       gene_descr <- summary$description,
       gene_desig <- summary$otherdesignations,
       gene_summary <- summary$summary
  )
}

expand_geneinfo <- function(gene_info_prep){
  gene_list <- lapply(gene_info_prep, summarize)
  return(gene_list)
}

#converting to dataframe 
convert_togenedf <- function(gene_list){
  gene_info_df <- data.frame(matrix(0, nrow = 1, ncol = 5))
  colnames(gene_info_df) <- c("name", "aliases", "description", "other designations", "summary")
  for (i in 1:length(gene_list))
  {
    gene_holder <- data.frame(gene_list[[i]])
    colnames(gene_holder) <- c("name", "aliases", "description", "other designations", "summary")
    
    gene_info_df<- rbind(gene_info_df, gene_holder)
  }
  
  return(gene_info_df)
}

#test --------------------------------------
#Object names are the same as those returned by the function called for simplicity.
#Can be different depending on user preferene.

#Define the search term
term <- "atrial fibrillation and cardiomyopathy"

#Gene search
gene_ids <- search_genes(term, retmax = 100)

#Obtain gene info
gene_info_prep <- obtain(gene_ids)

#Expand gene info
gene_list <- expand_geneinfo(gene_info_prep)

#Convert to dataframe and view
gene_info_df <- convert_togenedf(gene_list)
View(gene_info_df)

#LABELLING 

#function -------------------------------------------
# 1 = ion channel, 2 = sarcomeric, 3 = metabolic
afib_cm_label <- function(gene_info_df){
  for(i in 1:nrow(gene_info_df))
  {
    comparison <- paste(gene_info_df[i, ]$'description', gene_info_df[i, ]$'other designations', gene_info_df[i, ]$'summary', collapse = ",")
    gene_info_df$Classification[i] <- ifelse(grepl("channel|voltage", comparison), 1,
                                             ifelse (grepl("muscle|desmin|lamin|myosin|sarcomer|titin", comparison), 2,
                                                     ifelse (grepl("metabol|fatty|glycoly|citric acid|oxidati|mitochondr", comparison), 3,0)
                                             ))
  }
  return(gene_info_df)
}

#test ---------------------------------------------
gene_df_labelled <- afib_cm_label(gene_info_df)
View(gene_df_labelled)
