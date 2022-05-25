library(readxl)
library(data.table)

# Create data.Rdata file

tabela = "G:\\Meu Drive\\Doutorado\\tabela amostras.xlsm"

# load files
annotations_raw = setDT(read_excel(tabela, sheet="main"))
morphotypes_raw = setDT(read_excel(tabela, sheet="morfotipos"))


annotations_raw2 = annotations_raw[, c(3, 5:7)]
colnames(annotations_raw2) = c("cruzeiro", "operacao", "morfotipo", "Ninds")
morphotypes_raw = morphotypes_raw[, 1:4]


morphotypes_raw2 = fread("G:\\Meu Drive\\Doutorado\\morfotipos.csv")
morphotypes_raw2 = morphotypes_raw2[, c(2:6, 9)]


dados = na.omit(annotations_raw2)
dados = dados[, .(Ninds = sum(Ninds)), by = .(cruzeiro, operacao, morfotipo)]
dados = morphotypes_raw[dados, on = .(morfotipo)]
dados = morphotypes_raw2[dados, on = .(morfotipo)]


# filter
dados = dados[grupo != "ignorar"]
dados = dados[substr(operacao, 1, 2) != "BC"]

dados[grep("-", operacao), operacao := tstrsplit(operacao, "-")[[1]]]



# get environment data
distancia_rift = fread('draga_dist.csv')
dados = distancia_rift[dados, on = .(cruzeiro, operacao)]

envs = fread("draga_vars.csv")
dados = envs[dados, on = .(cruzeiro, operacao)]

# to factor
dados[, local := factor(local, c('North','South'))]
dados[grupo == "outros", grupo := "Others"]
dados[, grupo := factor(grupo, c('Porifera', 'Cnidaria', 'Annelida', 'Arthropoda', 'Echinodermata', 'Mollusca', 'Others'))]
dados[, cruzeiro := factor(cruzeiro, c("RGR1", "DY94"))]
dados[, operacao := factor(operacao)]

dados = dados[operacao != "RA1"]

dados[grupo == "Others", Class := Phylum]
dados[, Class := factor(Class, sort(unique(dados$Class))[c(8,11,1,20,13,17,21,19,23,14,7,2,16,9,12,22,18,3,10,15,4,5,6)])]


# Host table
hosp = annotations_raw[grepl("^[0-9]{1,}$", hospedeiro), c(3, 5:7, 11)]
colnames(hosp) = c("cruzeiro", "operacao", "morfotipo", "Ninds", "hospedeiro")
hosp = hosp[, .(Ninds = sum(Ninds)), by = .(cruzeiro, operacao, morfotipo, hospedeiro)]
hosp[, hospedeiro := as.integer(hospedeiro)]
hosp = morphotypes_raw[hosp, on = .(morfotipo)]

hosp[annotations_raw, hosp_morphotype := i.morfotipo, on = c(hospedeiro = "código")]

hosp = morphotypes_raw[hosp, on = c(morfotipo = "hosp_morphotype")]
hosp[, c("taxonRank", "hospedeiro") := NULL]
colnames(hosp)[1:7] = c("hosp_morfotipo", "hosp_taxon", "hosp_grupo", "morfotipo", "taxon", "taxonRank", "grupo")

hosp = hosp[grupo != "ignorar"]

save(dados, hosp, file = "data.Rdata")
