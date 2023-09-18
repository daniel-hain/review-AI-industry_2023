###########################################################################################
########################### Preamble
###########################################################################################

### Generic preamble
rm(list=ls())
set.seed(1337)

### Load packages  
library(tidyverse)
library(magrittr)

### Extra packages
# Biblis & NWs
library(bibliometrix)
library(tidygraph)
require(RNewsflow)
#NLP
library(tidytext)
library(topicmodels)
library(textstem)
library(ldatuning)

###########################################################################################
########################### Variable definitions
###########################################################################################

rm(list=ls())
# own Parameters
source("functions/00_parameters.R")
source("functions/functions_basic.R")
source("keys.R")

# institute and department
var_inst <- 'AI'
var_dept <- 'v1'

# openaikey
Sys.setenv(OPENAI_API_KEY = key_openai)

###########################################################################################
########################### Load & preprocessing articles
###########################################################################################

print('Starting: Loading Files')

files <- list.files(path = 'data', pattern = paste0('scopus_'), full.names = TRUE)

# Load bibliographic data
M <- convert2df(file = files, dbsource = "scopus", format = "csv") %>% 
  # Delete duplicates 
  distinct(UT, .keep_all = TRUE) 

# Filter 
M %<>% 
  # Filter number references
  mutate(CR_n = CR %>% str_count(';')) %>%
  # Abstract
  filter(AB != '') %>%
  filter(AB %>% str_length() >= 25) %>%
  # Number of cited references and citations
  mutate(TC_year = TC / (PY_max + 1 - PY)) %>%
  filter(TC_year >= 1)

#M_save <- M
#M <- M_save

# create label & Rownames
M %<>% rownames_to_column('XX') #%>% 
  #mutate(XX = paste(str_extract(XX, pattern = ".*\\d{4}"), str_sub(TI, 1,25)) %>% str_replace_all("[^[:alnum:]]", " ") %>% str_squish() %>% str_replace_all(" ", "_") %>% make.unique(sep='_'))
  
rownames(M) <- M$XX

# Save whole compilation
M %>% saveRDS(paste0('../temp/M_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))


# M <- read_rds(paste0('../temp/M_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))
M %>% biblioAnalysis(sep = ";") %>% saveRDS(paste0('../temp/M_res_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))

###########################################################################################
########################### Networks Bibliographic
###########################################################################################

print('Starting: Bibliographic Coupling network')

mat_bib <- M  %>% biblioNetwork(analysis = "coupling", network = "references", sep = ";", short = TRUE, shortlabel =  FALSE)
mat_bib %>% saveRDS(paste0('../temp/mat_bib__', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))
# mat_bib <- readRDS(paste0('../temp7mat_bib__', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))

g_bib <- mat_bib %>% igraph::graph_from_adjacency_matrix(mode = "undirected", weighted = TRUE, diag = FALSE) %>% 
  igraph::simplify() %>%
  as_tbl_graph(directed = FALSE) %N>% 
  left_join(M %>% select(XX, UT, PY, CR_n, TC_year), by = c("name" = "XX"))

## Restrict the network
g_bib <- g_bib %E>% 
  filter(weight >= cutof_edge_bib)

g_bib <- g_bib %N>%
  filter(!node_is_isolated()) %N>% # Could add:  | int_dept == TRUE
  mutate(dgr = centrality_degree(weights = weight)) %N>% 
  filter(dgr >= cutof_node_bib)

# Jaccard weighting
g_bib <- g_bib %E>% 
  mutate(weight_jac = weight / (.N()$CR_n[from] + .N()$CR_n[to] - weight) ) %E>%
  mutate(weight_jac = if_else(weight_jac > 1, 1, weight_jac) ) %N>%
  mutate(dgr_jac = centrality_degree(weights = weight_jac)) 

# # Further restrictions
# g_bib <- g_bib  %N>%
#   filter(percent_rank(dgr_jac) >= cutof_node_pct_bib) %E>% 
#   filter(percent_rank(weight_jac) >= cutof_edge_pct_bib) %N>%
#   filter(!node_is_isolated())

## Community Detection
g_bib <- g_bib %N>%
  mutate(com = group_louvain(weights = weight_jac)) %>%
  morph(to_split, com) %>% 
  mutate(dgr_int = centrality_degree(weights = weight_jac)) %N>%
  unmorph()

# Community size restriction
com_size_bib <- (g_bib %N>% as_tibble() %>% nrow()) * 0.05
g_bib %N>% as_tibble() %>% count(com, sort = TRUE)

g_bib <- g_bib %N>%
  group_by(com) %>%
    mutate(com_n = n()) %>%
  ungroup() %>%
  mutate(com = ifelse(com_n >= com_size_bib, com, NA) ) %>%
  mutate(com = ifelse(com <= com_max_bib, com, NA) ) %>%
  select(-com_n)  

# Delete nodes withou community
g_bib <- g_bib %N>%
  filter(!is.na(com))

# Update degree
g_bib <- g_bib %N>%
  mutate(dgr = centrality_degree(weights = weight),
         dgr_jac = centrality_degree(weights = weight_jac))


# Save the objects we need lateron
g_bib %>% saveRDS(paste0('../temp/g_bib_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))

## Merge with main data
M_bib <- M %>% select(UT, XX, PY) %>% inner_join(g_bib %N>% as_tibble() %>% select(UT, dgr, dgr_jac, com, dgr_int), by = 'UT') %>%
  distinct(UT, .keep_all = TRUE) 

M_bib %>% saveRDS(paste0('../temp/M_bib_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))


## Aggregated Network
g_bib_agg <- g_bib %N>%
  filter(!is.na(com)) %>%
  network_aggregate(by = "com", edge_attribute = "weight_jac", agg_FUN = sum)  %>%
  as.undirected(mode = "collapse", edge.attr.comb = "sum") %>%
  as_tbl_graph(directed = FALSE) %N>%
  select(-name) %>%
  mutate(id = 1:n()) %E>%
  rename(weight = agg.weight_jac) %>%
  select(from, to, weight)

## Weight edges
# g_bib_agg <- g_bib_agg %E>%
#   rename(weight_count = weight) %>%
#   mutate(weight = weight_count / (.N()$N[from] * .N()$N[to]) ) %>%
#   mutate(weight = (weight * 100) %>% round(4)) %N>%
#   mutate(dgr = centrality_degree(weights = weight))

# Save the objects we need lateron
g_bib_agg %>% saveRDS(paste0('../temp/g_bib_agg_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))

# Delete all we dont need
rm(mat_bib, g_bib, com_size_bib, cutof_edge_bib, cutof_node_bib, g_bib_agg)




###########################################################################################
########################### Collaboration NW
###########################################################################################

pub_inst <- M %>% as_tibble() %>% metaTagExtraction(Field = "AU_UN") %>% 
  select(XX, UT, PY, AU_UN) %>% 
  # Sepperate AU_UN for 2_m edgelist
  separate_rows(AU_UN, sep = ';') %>%
  # filter
  drop_na(AU_UN) %>%
  filter(!(AU_UN %in% c('', ' ', 'NA', 'NOTREPORTED', 'NOTDECLARED'))) %>%
  # Only 1 link per paper, independent of author number
  distinct(UT, AU_UN, .keep_all = TRUE) 

el_inst <- pub_inst %>% 
  left_join(pub_inst %>% select(UT, AU_UN), by = 'UT') %>% 
  rename(from = AU_UN.x, to = AU_UN.y) %>%
  filter(from != to) %>%
  group_by(from, to) %>%
  summarise(weight = n()) %>%
  ungroup()

# Save it 
pub_inst %>% saveRDS(paste0('../temp/pub_inst_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))
el_inst %>% saveRDS(paste0('../temp/el_inst_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))

# clean 
# rm(pub_inst, el_inst)

###########################################################################################
########################### Network Cocitation 
###########################################################################################

print('Starting: Co-Citation Network')

mat_cit <- M %>%
  semi_join(M_bib, by = 'UT') %>%
  as.data.frame() %>% 
  biblioNetwork(analysis = "co-citation", network = "references", sep = ";", short = TRUE, shortlabel = FALSE)

mat_cit %>% saveRDS(paste0('../temp/mat_cit_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))
# mat_cit <- readRDS(paste0('../temp/mat_cit_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))

g_cit <- mat_cit %>% igraph::graph_from_adjacency_matrix(mode = "undirected", weighted = TRUE, diag = FALSE) %>% 
  igraph::simplify() %>%
  as_tbl_graph(directed = FALSE) # %N>% left_join(M %>% select(XX, SR, PY, TC, J9), by = c("name" = "XX")) %>% mutate(id = 1:n())

# Restrict the network
g_cit <- g_cit %E>% 
  filter(weight >= cutof_edge_cit) %N>%
  filter(!node_is_isolated())

g_cit <- g_cit %N>%
  mutate(dgr = centrality_degree(weights = weight)) %N>%
  filter(dgr >= cutof_node_cit) 

# Further restrictions
g_cit <- g_cit %N>% 
  filter(percent_rank(dgr) >= cutof_node_pct_cit) %E>%
  filter(percent_rank(weight) >= cutof_edge_pct_cit) %N>%
  filter(!node_is_isolated())

## Community Detection
g_cit <- g_cit %N>%
  mutate(com = group_louvain(weights = weight)) %N>%
  morph(to_split, com) %>% 
  mutate(dgr_int = centrality_degree(weights = weight)) %>%
  unmorph()

g_cit %N>% as_tibble() %>% count(com)

# community detection
com_size_cit <- (g_cit %N>% as_tibble() %>% nrow()) * 0.05

# Community size restriction
g_cit <- g_cit %N>%
  group_by(com) %>%
  mutate(com_n = n()) %>%
  ungroup() %>%
  mutate(com = ifelse(com_n >= com_size_cit, com, NA) ) %>%
  mutate(com = ifelse(com <= com_max_cit, com, NA) ) %>%
  select(-com_n)  

# Delete nodes withou community
g_cit <- g_cit %N>%
  filter(!is.na(com))

# Update degree
g_cit <- g_cit %N>%
  mutate(dgr = centrality_degree(weights = weight))

# Save the objects we need lateron
g_cit %>% saveRDS(paste0('../temp/g_cit_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))

# generate citation report
C_nw <- g_cit %N>% as_tibble() 
C_nw %>%  saveRDS(paste0('../temp/C_nw_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))

## A ggregated Network
require(RNewsflow)
g_cit_agg <- g_cit %>%
  network_aggregate(by = "com", edge_attribute = "weight", agg_FUN = sum)  %>%
  as.undirected(mode = "collapse", edge.attr.comb = "sum") %>%
  as_tbl_graph(directed = FALSE) %N>%
  select(-name) %>%
  mutate(id = 1:n()) %E>%
  rename(weight = agg.weight) %>%
  select(from, to, weight)

g_cit_agg %>% saveRDS(paste0('../temp/g_cit_agg_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))

rm(mat_cit, g_cit, g_cit_agg)

###########################################################################################
########################### 2 mode network 
###########################################################################################

print('Starting: 2 Mode Network')

M <- read_rds(paste0('../temp/M_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))
M_bib <- read_rds(paste0('../temp/M_bib_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))
C_nw <- read_rds(paste0('../temp/C_nw_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))

#rownames(M) <- M %>% pull(UT)

m_2m <- M %>% 
  semi_join(M_bib, by = 'UT') %>%
  as.data.frame() %>% cocMatrix(Field = "CR", sep = ";", short = TRUE)

g_2m <- m_2m %>% igraph::graph_from_incidence_matrix(directed = TRUE, mode = 'out', multiple = FALSE) %>% 
  igraph::simplify() 

el_2m <- g_2m %>%
  get.edgelist() %>%
  as_tibble() %>%
  rename(from = V1,
         to = V2)

el_2m %<>%
  left_join(M_bib %>% select(XX, com, PY), by = c('from' = 'XX')) %>%
  rename(com_bib = com) %>%
  left_join(C_nw %>% select(name, com), by = c('to' = 'name')) %>%
  rename(com_cit = com) %>% 
  drop_na(PY, com_bib, com_cit)

# save
el_2m %>% saveRDS(paste0('../temp/el_2m_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))

rm(m_2m, g_2m, el_2m, C_nw)


###########################################################################################
########################### NLP - Research Areas
########################################################################################### 

# Extract all for externbal BERTopic Modelling
df_text <- M %>% 
  as_tibble() %>%
  select(UT, PY, TI, AB) %>%
  left_join(M_bib %>% select(UT, com, dgr_int), by = 'UT') %>%
  mutate(text = paste(TI, AB, sep = '. ')  %>% 
           str_to_lower() %>% 
           str_remove_all("©.*") %>%          
           str_remove_all("/(&trade;|&reg;|&copy;|&#8482;|&#174;|&#169;)/.*") %>%
           str_remove_all('(&elsevier;|&springer;|&rights reserved)/.*') %>%
           str_squish())  %>%
  select(UT, PY, text, com, dgr_int) 

df_text %>% write_csv('data/data_text_all.csv')

df_text <- read_csv('data/data_text_all.csv')

library(openai)

M_bib %>% count(com)

top_n_bib = 10

promt_intro = 
"You are a knowledgable and helpful researcher in social science. I want you to summarize the research in the following documents."

promt_context = "All the documents are more or less strongly related to and cover diverent aspects of research on the impact of 
artificial intelligence on industry dynbamics and innovation.
Therefore, they all should be interpreted in relation to this overal theme."

### 

promt_bib_intro = 
  "I will provide you text with titles plus abstacts of scientific journal article publications. They are representative articles for broader research themes to be identified. 
They are supposed to be similar in terms of theoretical foundations, literature they relate to, context or topic of research."

promt_bib_instruction = 
  "Your task is to summarize the topic overal research based on the provided article text by a short label of 2-5 words, plus a short description of 3-5 sentences. 
The label optimally relates to relvevant scientific concepts in the field of study, and maybe the context it is studied.
Only create a single label and description summarizing the overarching research theme in the documents.
The description should be in light of the overal context, brief, focussed, clear, and avoid redundancies. This summary should highlight the commonality of the documents. 
It should indicate the main theoretical theme, research framework applied, context, potential contributions and implications for theory, policy, and professionals."
# You response should have the following format: <label> | <description>" 

promt_docs_bib = paste("I will now provide the", top_n_bib, "documents. Every document starts with an '-', and ends with a linebreak.", sep = " ")

promt_bib_doc <- df_text %>% 
  filter(com == 6) %>%
  slice_max(order_by = dgr_int, n = top_n_bib) %>%
  pull(text) %>% paste('-', ., sep = ' ', collapse = ' \n ')

promt_bib <- paste(promt_bib_intro, promt_context, promt_bib_instruction, promt_docs_bib, promt_bib_doc, sep = ' \n \n ') # 
promt_bib

desc_bib <- create_chat_completion(
  model = "gpt-3.5-turbo",
  messages = list(
    list(
      "role" = "system",
      "content" = promt_intro
    ),
    list(
      "role" = "user",
      "content" = promt_bib
    )
  )
)



###########################################################################################
########################### NLP - Research Areas
########################################################################################### 

C_nw <- read_rds(paste0('../temp/C_nw_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))

C_nw %>% count(com)

top_n_cit = 50

promt_cit_intro = 
  "I will provide you are titles of scientific journal article publications in the following format: Author, Publication title (Publication year). 
They are representative articles for broader stream of literature to be identified. 
the articles are references which are often cited together (co-citation), therefore often represent seminal articles of a research field."

promt_cit_instruction = 
  "Your task is to summarize the topic of overal research based on the provided article text by a short label of 2-5 words, plus a short description of 3-5 sentences. 
The label optimally relates to relvevant scientific concepts in the field of study, and maybe the context it is studied.
Only create a single label and description summarizing the overarching research theme in the documents.
The description should be in light of the overal context, brief, focussed, clear, and avoid redundancies. This summary should highlight the commonality of the documents. 
It should indicate the main theoretical theme, research framework applied. and main arguments made."
# You response should have the following format: <label> | <description>" 

promt_docs_cit = paste("I will now provide the", top_n_cit, "documents. Every document starts with an '-', and ends with a linebreak.", sep = " ")

promt_cit_doc <- C_nw %>% 
  filter(com == 6) %>%
  slice_max(order_by = dgr_int, n = top_n_cit) %>%
  pull(name) %>% paste('-', ., sep = ' ', collapse = ' \n ')

promt_cit <- paste(promt_cit_intro, promt_context, promt_cit_instruction, promt_docs_cit, promt_cit_doc, sep = ' \n \n ') # 
promt_cit

desc_cit <- create_chat_completion(
  model = "gpt-3.5-turbo",
  messages = list(
    list(
      "role" = "system",
      "content" = promt_intro
    ),
    list(
      "role" = "user",
      "content" = promt_bib
    )
  )
)


###########################################################################################
########################### NLP Topics
########################################################################################### 

df_top <- read_csv('data/documents_topic.csv')[-1] %>% 
  as_tibble() %>% 
  pivot_longer(where(is.numeric), names_to = 'com', values_to = 'weight') %>%
  mutate(com = as.numeric(com)) %>%
  filter(weight >= 0.05)

M_top <- M %>% select(XX, UT, PY, TI, SO, TC_year) %>%
  inner_join(df_top, by = 'UT')

M_top %>% write_rds(paste0('../temp/M_top_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))

df_text_top <- M %>% 
  as_tibble() %>%
  select(UT, PY, TI, AB) %>%
  inner_join(M_top %>% select(UT, com, weight), by = 'UT') %>%
  mutate(text = paste(TI, AB, sep = '. ')  %>% 
           str_to_lower() %>% 
           str_remove_all("©.*") %>%          
           str_remove_all("/(&trade;|&reg;|&copy;|&#8482;|&#174;|&#169;)/.*") %>%
           str_remove_all('(&elsevier;|&springer;|&rights reserved)/.*') %>%
           str_squish())  %>%
  select(UT, PY, text, com, weight) 

M_top %>% count(com)

top_n_top = 10

promt_top_intro = 
  "I will provide you text with titles plus abstacts of scientific journal article publications. They are representative articles for broader research topic to be identified. 
They are supposed to be similar in terms of theoretical foundations, literature they relate to, context or topic of research."

promt_top_instruction = 
  "Your task is to summarize the topic overal research based on the provided article text by a short label of 2-5 words, plus a short description of 3-5 sentences. 
The label optimally relates to relvevant scientific concepts in the field of study, and maybe the context it is studied.
Only create a single label and description summarizing the overarching research theme in the documents.
The description should be in light of the overal context, brief, focussed, clear, and avoid redundancies. This summary should highlight the commonality of the documents. 
It should indicate the main theoretical theme, research framework applied, context, potential contributions and implications for theory, policy, and professionals.
You response should have the following format: <label> | <description>" 
#

promt_docs_top = paste("I will now provide the", top_n_top, "documents. Every document starts with an '-', and ends with a linebreak.", sep = " ")

promt_top_doc <- df_text_top %>% 
  filter(com == 3) %>%
  slice_max(order_by = weight, n = top_n_top, with_ties = FALSE) %>%
  pull(text) %>% paste('-', ., sep = ' ', collapse = ' \n ')

promt_top <- paste(promt_top_intro, promt_context, promt_top_instruction, promt_docs_top, promt_top_doc, sep = ' \n \n ') # 
promt_top

# # TODO: TRYOUT THIS FUINCTION FOR AUTOMATIOZATION
# ALSO; SOLVE LEAKED API KEY ISSUE
# promt_system_top <- promt_intro
# promt_context_top <- paste(promt_top_intro, promt_context, promt_top_instruction, promt_docs_top, sep = ' \n \n ') 
# 
# res_top <- create_gpt_summary(df_text = df_text_top, 
#                               promt_system = promt_system_top, 
#                               promt_context = promt_context_top, 
#                               select_type = 'TP', 
#                               top_n = 10,
#                               model_version = "gpt-4")

###########################################################################################
########################### Local citations
########################################################################################### 

print('Starting: Further Analysis')

CR <- M %>% citations(sep = ";")
CR %>% saveRDS(paste0('../temp/CR_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))

#CRL <- M %>% localCitations(sep = ";") # For some reason takes forever...
#CRL %>% saveRDS(paste0('../temp/CRL_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))

rm(CR)

###########################################################################################
############################ Threefield Plot
########################################################################################### 

M_threefield <- M %>% as.data.frame() %>% threeFieldsPlot(fields = c("AU_UN", "DE", "CR_SO"), n = c(20, 10, 10))
M_threefield %>% saveRDS(paste0('../temp/threefield_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))
rm(M_threefield)
