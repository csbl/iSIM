######## Read and save the iSIM model from TSV files ########
# This section only needs to be run once  
# Download packages necessary for this analysis
source("http://bioconductor.org/biocLite.R")
# sybil: necessary to run FVA and FBA on metabolic models
biocLite("sybil")
# glpkAPI: necessary to run GLPK-based optimization in the sybil package
biocLite("glpkAPI")
# dplyr: useful for manipulating dataframe
biocLite("dplyr")
# reshape2: useful for reshaping data frames and matrices
biocLite("reshape2")

######## Load packages ########
library(reshape2)
library(dplyr)
library(sybil)
library(glpkAPI)

# Read in the sybil version of iSIM (takes a minute or two)
# read in model in tabular form
iSIM = readTSVmod(prefix = "iSIM", "tsv", balanceReact = F, def_bnd = 1000)


model.metabolite.info = data_frame(
  met_id = iSIM@met_id,
  met_name = iSIM@met_name)
# Check out metabolites in the model
View(model.metabolite.info)

model.reaction.info = data_frame(
  react_id = iSIM@react_id,
  react_name = iSIM@react_name,
  react_objective = iSIM@obj_coef)
# Check out reactions in the model
View(model.reaction.info)

# Save the sybil model as a file to avoid re-loading model.
saveRDS(iSIM,"iSIM.rda")

#!!!!!!!!!!!!!!!!!!! Start here if you have already run the previous section !!!!!!!!!!!!!!!!!!!#

######## Load iSIM and necessary packages ########

# Load packages
library(dplyr)
library(reshape2)
library(sybil)
library(glpkAPI)

# Load the sybil model (saved from readTSVmod)
iSIM = readRDS("iSIM.rda")

######## Get reaction info and make new reaction formula strings ########
# This section may be difficult to understand without knowledge of the dplyr package
# "%>%" means "then...", see dplyr package documentation

# Reaction formulas are not represented in a sybil model except in the S-matrix
iSIM_s_matrix = as.matrix(iSIM@S)
rownames(iSIM_s_matrix) = iSIM@met_id
colnames(iSIM_s_matrix) = iSIM@react_id

iSIM_rxn_decomposition = iSIM_s_matrix %>% 
  # transform the S matrix ( met x rxn matrix with coef as value ) 
  # to the long form of a reaction equation (one reaction represented as multiple rows)
  melt %>% # then...
  # rename default variables output by melt
  transmute(met_id = as.character(Var1), rxn_id = as.character(Var2), coef = as.numeric(value)) %>% # then...
  # remove all of the zero value coefficients from the S matrix
  filter(coef != 0) %>% # then...
  # extract compartment and (compartmentless) metabolic compound ids
  # from the met_id which is in the form: "m[compound_id_number][compartment_abbreviation]"
  mutate(met_compartment = gsub("m[0-9]+","",met_id), 
         met_compound = gsub("[a-z]$","",met_id)) %>% as.tbl %>%
  # pair met_id with the cpd_name from the original model
  left_join(data_frame(met_id = iSIM@met_id, cpd_name = iSIM@met_name)) %>%
  # new compound name includes compartment of metabolite
  mutate(met_name = paste0(cpd_name, "[",met_compartment,"]"))

# create the S matrix by replacing row and column names with metabolite (full name with compartment) and reaction names
iSIM@S %>% as.matrix %>% 
  # extract only the compartment ID from the met_id
  data.frame(row.names = paste0(iSIM@met_name, "[",gsub(".*([a-z])$","\\1",iSIM@met_id),"]")) %>% 
  # change the names of the columns
  setNames(iSIM@react_name) %>% 
  # save the table
  write.table("iSIM_smatrix_name.txt",quote = F,sep = "\t",row.names = T,col.names = T)

# create the S matrix by replacing row and column names with metabolite IDs and reaction IDs
iSIM@S %>% as.matrix %>% data.frame(row.names = iSIM@met_id) %>% setNames(iSIM@react_id) %>% 
  write.table("iSIM_smatrix_id.txt",quote = F,sep = "\t",row.names = T,col.names = T)


# display the first few rows:
iSIM_rxn_decomposition  %>% head

# create new rxn formulas given the reaction matrix decomposition
iSIM_rxn_info = iSIM_rxn_decomposition %>% 
  # annotate metabolite names from iSIM sybil model (by joining tables)
  left_join(data_frame(met_id = iSIM@met_id,met_abbrev = iSIM@met_name),by="met_id") %>%
  mutate(
    #cpd_base = ef_met_compound(met_id),
    # metabolites will be displayed by their name + compartment
    met_base = paste0(cpd_name,"[",met_compartment,"]"),
    # define metabolites as either substrates (left side of reaction) or products (right side of reaction)
    met_substrates = coef<0,met_products = coef>0,
    # append stoichiometric coefficients to metabolites when coef != 1
    coef_base = ifelse(abs(coef) == 1,"",paste0(as.character(abs(coef))," ")),
    met_coef_base = paste0(coef_base,met_base)) %>% 
  # group by each reaction then paste together each side of the reaction equation (substrates, products)
  group_by(rxn_id) %>% summarize(
    base_substrates = paste(met_coef_base[met_substrates],collapse = " + "),
    base_products = paste(met_coef_base[met_products],collapse = " + ")) %>% ungroup %>%
  # annotate reactions by their lower and upper bounds to determine directionality of the arrow
  inner_join(data_frame(rxn_id = iSIM@react_id,
                        rxn_name = iSIM@react_name,
                        lb = iSIM@lowbnd,
                        ub = iSIM@uppbnd) %>%
               mutate(base_arrow = ifelse(lb < 0," <==> "," --> ")), by="rxn_id") %>%
  # paste together substrates and products by the reaction arrow to generate fresh reaction formulas
  mutate(base_substrates = ifelse(is.na(base_substrates)," ",base_substrates),
         base_products = ifelse(is.na(base_products)," ",base_products),
         rxn_formula = paste0(base_substrates,base_arrow,base_products)) %>% 
  # remove temporary variables (start with base_)
  select(-base_substrates,-base_products,-base_arrow)

# display the first few rows:
iSIM_rxn_info %>% head 

#### Run flux balance analysis on the default model ####

# perform FBA on the model as loaded (same parameters as cobra model I had sent in November 2014)
# the Default model's objective is to maximize ATP hydrolysis:
iSIM_rxn_info %>% 
  filter(rxn_id %in% iSIM@react_id[iSIM@obj_coef > 0]) %>% select(rxn_id,rxn_formula)
# In Ratcon1, RCR11017: ATP[c] + H20[c] --> ADP[c] + Pi[c]
# This reaction simulates the consumption of ATP for energy-dependent
# processes such as active transporter or actin-myosin contraction


# FIGURE: colored stoichiometric matrix
library(ggplot2)
plot.smatrix = iSIM_rxn_decomposition %>% 
  left_join(iSIM_rxn_info) %>% 
  mutate(met_name = factor(met_name, levels = rev(unique(met_name)), ordered = T)) %>% 
  arrange(grepl("E",rxn_id), rxn_id) %>% mutate(rxn_name = factor(rxn_name, levels = unique(rxn_name), ordered = T)) %>% 
  # group_by(rxn_id) %>% mutate(value = coef / mean(abs(coef))) %>% ungroup %>% 
  ggplot(aes(x = rxn_name, y = met_name, fill = sign(coef), label = coef)) + 
  geom_tile() + 
  geom_text(size = 3, color = "#FFFFFF", fontface = "bold") + 
  scale_fill_gradient2(low = "#4F81BD", mid = "#FFFFFF", high = "#C0504D") + 
  theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") + 
  xlab(NULL) + ylab(NULL) 


plot.smatrix
ggsave(filename = "iSIM_smatrix.pdf",plot.smatrix,width = 2.6, height = 3)

# function to run FVA with inputs of a required percentage, model, and table to add results to
ef_tbl_fva = function(x.pct,x.model,x.tbl) {
  x.fva = fluxVar(x.model,percentage = x.pct,fixObjVal = T)
  x.result = lp_obj(x.fva) %>% round(9) 
  x.n = length(x.result)
  data_frame(rxn_id = x.fva@react@react_id,
             # fluxVar returns one vector with first half containing lb and second half containing ub
             fva_lb = x.result[1:(x.n/2)],
             fva_ub = x.result[(x.n/2+1):(x.n)]) %>% left_join(x.tbl) %>%
    mutate(fva_pct = x.pct,
           fva_req = fva_lb > 1e-9 | fva_ub < -1e-9,
           fva_on = abs(fva_lb) > 1e-9 | abs(fva_ub) > 1e-9) %>% as.tbl
}


# Run FVA while requiring 0-100% of the maximum ATP yield (increments of 5%)
fva.pct.result = lapply(c(0:20)*5,ef_tbl_fva,iSIM,iSIM_rxn_info) %>% bind_rows
# Display FVA results table
fva.pct.result %>% arrange(rxn_id)
# Write FVA results to a text file
fva.pct.result %>% write.table("iSIM_fva_result_percentage.txt", sep = "\t", quote = F, row.names = F)

# Prepare FVA results for plotting so that each lb or ub value 
# is represent in a different row (see reshape2::melt documentation)
fva.pct.data = bind_rows(list(
  fva.pct.result %>% select(rxn_id,rxn_name,fva_pct,fva_req,fva_on,lb = fva_lb,ub = fva_ub) %>%
    mutate(ub = ifelse(fva_on,ub,0.00001),smaller = ifelse(abs(lb)<abs(ub),"lb","ub")) %>% 
    melt(c("rxn_id","rxn_name","fva_pct","smaller","fva_req","fva_on")) %>% as.tbl %>% 
    mutate(lrg_bnd = ifelse(smaller == variable,"s","c"),bnd = "fva")))

# Plot FVA results of lower and upper bound ranges showing the tradeoff
# between the "glycolysis" and "respiration" reaction pathways
fva.pct.plot = fva.pct.data %>% filter(rxn_id %in% c("R1","R2")) %>% 
  mutate(color = ifelse(fva_req, "#C0504D", ifelse(fva_on, "#7F7F7F", "#4F81BD"))) %>%
  mutate(shape = ifelse(smaller == variable, 15, 18)) %>% 
  ggplot(aes(x = value, y = fva_pct, group = factor(fva_pct), color = color)) + 
  geom_line() + geom_vline(linetype = "solid",size = 1,xintercept =0) +
  # geom_point(size = 5,aes(shape = factor(lrg_bnd, levels= c("s", "c", "t")))) +
  geom_point(size = 3,aes(shape = shape)) +
  scale_shape_identity() + #(values = c("circle","square")) + 
  scale_color_identity() + 
  scale_x_continuous(breaks = c(0:10) / 5) + 
  scale_y_continuous(breaks = c(0:20)*5) + 
  theme_minimal(base_size = 16) + 
  facet_wrap("rxn_name") +
  theme(legend.position = "none") + xlab("Units of Flux through Reaction") + ylab("Required Flux through Obj Fun")

fva.pct.plot
ggsave(filename = "iSIM_fva_percentage.pdf",fva.pct.plot,width = 6, height = 3.5)

# Do the same as above except for increments of 1 ATP per step up to 32 (max)
fva.increment.result = lapply(100 * c(0:16) / 16,ef_tbl_fva,iSIM,iSIM_rxn_info) %>% bind_rows
# Write FVA results to a text file
fva.increment.result %>% write.table("iSIM_fva_result_increment.txt", sep = "\t", quote = F, row.names = F)

fva.increment.data = bind_rows(list(
  fva.increment.result %>% select(rxn_id,rxn_name,fva_pct,fva_req,fva_on,lb = fva_lb,ub = fva_ub) %>%
    mutate(ub = ifelse(fva_on,ub,0.00001),smaller = ifelse(abs(lb)<abs(ub),"lb","ub")) %>% 
    melt(c("rxn_id","rxn_name","fva_pct","smaller","fva_req","fva_on")) %>% as.tbl %>% 
    mutate(lrg_bnd = ifelse(smaller == variable,"s","c"),bnd = "fva")))

fva.increment.plot = fva.increment.data %>% filter(rxn_id %in% c("R2","E2")) %>% 
  mutate(color = ifelse(fva_req, "#C0504D", ifelse(fva_on, "#7F7F7F", "#4F81BD"))) %>%
  mutate(shape = ifelse(smaller == variable, 15, 18)) %>% 
  ggplot(aes(x = value, y = .32 * fva_pct, group = factor(fva_pct), color = "#C0504D")) + 
  geom_line(size = 1) + geom_vline(linetype = "solid",size = 1, xintercept = 0) +
  geom_point(size = 3,aes(shape = shape)) +
  scale_shape_identity() + 
  scale_color_identity() + 
  #scale_x_continuous(breaks = c(0:10) / 5) + 
  #scale_y_continuous(breaks = 0) + # breaks = 2 * c(0:16)) + 
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1)) + 
  facet_wrap("rxn_name", scales = "free_x") +
  theme(title = element_text(size = 12, face = "bold"), axis.title = element_text(size = 12, face = "bold"), strip.text = element_text(face = "bold")) +
  xlab("Flux range") + ylab("ATP demand requirement") + 
  theme(panel.background = element_rect(fill = "white"))

fva.increment.plot
ggsave(filename = "iSIM_fva_increment.pdf",fva.increment.plot,width = 5, height = 3)

# Plot the results for every reaction in the network
fva.increment.plot.all = fva.increment.data %>% #filter(rxn_id %in% c("E1","E3")) %>% #,!grepl("E", rxn_id)) %>%
  mutate(color = ifelse(fva_req, "#C0504D", ifelse(fva_on, "#7F7F7F", "#4F81BD"))) %>%
  mutate(shape = ifelse(smaller == variable, 15, 18)) %>% 
  # filter(fva_pct < 20) %>%
  ggplot(aes(x = value, y = .32 * fva_pct, group = factor(fva_pct), color = color)) + 
  geom_line() + geom_vline(linetype = "solid",size = 1, xintercept = 0) +
  # geom_point(size = 5,aes(shape = factor(lrg_bnd, levels= c("s", "c", "t")))) +
  geom_point(size = 3,aes(shape = shape)) +
  scale_shape_identity() + #(values = c("circle","square")) + 
  scale_color_identity() + 
  # scale_x_continuous(breaks = -(12:0)/.5) + 
  scale_y_continuous(breaks = 2 * c(0:16)) + 
  theme_minimal(base_size = 12) + 
  facet_wrap("rxn_name", scales = "free_x") +
  theme(legend.position = "none") + xlab(NULL) + ylab(NULL)

fva.increment.plot.all
ggsave(filename = "iSIM_fva_increment_all.pdf",fva.increment.plot.all,width = 10, height = 6)


#### Perform single gene deletions for all gene-associated reactions ####
# iSIM includes representative genes for several simplified reactions
# Note that exchange reactions do not typically have gene associations
# as we used here to keep things minimal.

gene.ko = iSIM %>% geneDeletion()
# delete every gene and record new FBA value
gene.ko.fba = gene.ko %>% fluxdels %>% setNames(gene.ko@dels) %>% sapply(function(x) {
  (iSIM %>% rmReact(x) %>% optimizeProb("fba",retOptSol = F))$obj
})

tbl.ko = data_frame(gene_id = names(gene.ko.fba),gene_ko_atp = gene.ko.fba) %>% left_join(
  data_frame(gene_id = gene.ko@dels[,1],gene_rxn_ko = gene.ko@fluxdels %>% 
               lapply(unlist) %>% lapply(unique) %>% lapply(paste0,collapse = ";") %>% sapply(c))) %>% 
  left_join(gene.ko %>% fluxdels %>% setNames(gene.ko@dels) %>% melt %>% transmute(rxn_id = as.character(value), gene_id = L1)) %>%
  left_join(iSIM_rxn_info)

View(tbl.ko)

# write table which includes FBA values for knockout with reaction information
tbl.ko %>% write.table("iSIM_gene_knockout_screen", sep = "\t", quote = F, row.names = F)

# perform double knockouts
gene.ko2 = iSIM %>% geneDeletion(combinations = 2)
# run the FBA simulation for each pair of gene knockouts
gene.ko2.fba = gene.ko2 %>% fluxdels %>% setNames(gene.ko2@dels %>% apply(1,paste0,collapse = "_")) %>% sapply(function(x) {
  (iSIM %>% rmReact(x) %>% optimizeProb("fba",retOptSol = F))$obj
})

# summarize the genes and reactions removed in a table
gene.ko2.rxns = gene.ko2@fluxdels %>% 
  setNames(names(gene.ko2.fba)) %>% melt %>% 
  transmute(rxn_id = as.character(value), genes = L1) %>%
  left_join(iSIM_rxn_info) %>% group_by(genes) %>% 
  summarize(rxn1 = first(rxn_id), rxn2 = last(rxn_id), 
            name1 = first(rxn_name), name2 = last(rxn_name), 
            rxns = paste0(rxn_id, collapse = "_"), 
            names = paste0(rxn_name, collapse = " / ")) %>% ungroup


tbl.ko2 = data_frame(genes = names(gene.ko2.fba),atp = gene.ko2.fba) %>% 
  left_join(gene.ko2.rxns) %>%
  # insert the ATP values for individual gene KOs
  left_join(tbl.ko %>% select(rxn1 = rxn_id, atp1 = gene_ko_atp)) %>%
  left_join(tbl.ko %>% select(rxn2 = rxn_id, atp2 = gene_ko_atp)) %>%
  mutate(atp12 = atp)
tbl.ko2 %>% filter(atp12 < atp1, atp12 < atp2) %>% View
