library(readxl)
library(data.table)
library(vegan)
library(ggplot2)
# library(svglite)

load("data.Rdata")


# Plot - abundance vs class
tmp = dados[, .(Ninds = sum(Ninds), S = length(unique(morfotipo))), by = .(grupo, Class)]

g <- ggplot(dados, aes(x=Class, y=Ninds)) + geom_col() +
  geom_text(aes(label = S), tmp, vjust = -0.3, size = 9/ggplot2:::.pt) +
  facet_grid(~ grupo, 
             scales = "free_x", # Let the x axis vary across facets.
             space = "free_x",  # Let the width of facets vary and force all bars to have the same width.
             switch = "x") +
  scale_y_continuous(expand = expansion(c(0, 0.05))) +
  ylab("Abundance") +
  theme_classic(base_size = 8, base_line_size = 0.25) +
  theme(strip.placement = "outside",                      # Place facet labels outside x axis labels.
        strip.background = element_blank(),  # Make facet label background white.
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1))

g

ggsave("abundance.png", g, height = 12, width = 16, units = 'cm', dpi = 600)




# most abundant morphotypes - used in the results
x=dados[, .(Ninds=sum(Ninds)), by=.(morfotipo)]
setorder(x, -Ninds)
head(x, 10)



##########################
# prepare data 
# remove sites too far away from the rift
dados2 = dados[!(cruzeiro == 'RGR1' & operacao %in% c('D02', 'D03', 'D04', 'D11', 'D12'))]
dados2 = dados2[!(cruzeiro == 'DY94' & operacao == 'HY39')]


envs = dados2[, c(.(Ninds = sum(Ninds)), lapply(.SD, first)), by = .(cruzeiro, operacao), .SDcols = c('dist_rift', 'crusts', 'shape_leng', 'local', 'depth', 'slope', 'mean_dir')]
setorder(envs, cruzeiro, operacao)

mc = dcast(dados2, cruzeiro + operacao ~ morfotipo, value.var = 'Ninds', fun = sum)
setorder(mc, cruzeiro, operacao)
mc[, c('cruzeiro', 'operacao') := NULL]

envs$S <- specnumber(mc)
envs$H <- diversity(mc, index = "shannon")
envs$J <- envs$H / log(envs$S)
envs[, op := ifelse(cruzeiro == 'RGR1', 'D1', ifelse(substr(operacao , 1, 1) == "D", "D2", "HY"))]


# correlation between abundance and shape_leng
cor.test(envs$shape_leng, envs$Ninds)
cor.test(envs$shape_leng, envs$S)
# plot(dados2$shape_leng, dados2$Ninds)


anova(lm(Ninds ~ op, envs))
anova(lm(S ~ op, envs))

## permanova
mc2 = log1p(mc)

mod <- betadisper(vegdist(mc2), envs$op)
permutest(mod, permutations = 9999)

fit <- adonis2(mc2 ~ op + shape_leng, envs)
fit


### difference in community based in op, try with only dredges

###########################
dados2 = dados2[substr(operacao , 1, 1) == "D"] # remove HyBIS and leave only dredges


envs = dados2[, c(.(Ninds = sum(Ninds)), lapply(.SD, first)), by = .(cruzeiro, operacao), .SDcols = c('dist_rift', 'crusts', 'shape_leng', 'local', 'depth', 'slope', 'mean_dir')]
setorder(envs, cruzeiro, operacao)

mc = dcast(dados2, cruzeiro + operacao ~ morfotipo, value.var = 'Ninds', fun = sum)
setorder(mc, cruzeiro, operacao)
mc[, c('cruzeiro', 'operacao') := NULL]



# correlation between abundance and shape_leng
cor.test(envs$shape_leng, envs$Ninds)
cor.test(envs$shape_leng, envs$S)
# plot(dados2$shape_leng, dados2$Ninds)


anova(lm(Ninds ~ cruzeiro, envs))
anova(lm(S ~ cruzeiro, envs))

## permanova
mc2 = decostand(mc, "hellinger") #log1p(mc)
mcd = dist(mc2) # vegdist(mc2)

mod <- betadisper(mcd, envs$cruzeiro)
permutest(mod, permutations = 9999)

  

### now the permanova test is not significant
### Continue with further analysis



# check collinearity
library(usdm)
M <- vifcor(envs[,.(dist_rift, crusts, depth, slope, mean_dir)])
M



# permanova
mod <- betadisper(mcd, envs$local)
permutest(mod, permutations = 999)

fit <- adonis2(mcd ~ local + dist_rift + crusts + depth + slope + mean_dir, envs, by = "terms", permutations = 999)
fit


# NMDS
set.seed(12345)
NMDS <- metaMDS(mcd, k=2, try = 1000, trymax = 10000)
depth = envs$depth * -1
fit_env = envfit(NMDS ~ depth)
localf = as.numeric(envs$local)


library(fplot)

# or pdf_fit
png_fit("nmds.png", pt = 8, width = "14cm", height = "8.5cm", res= 1000)
par(mar=c(2.5,2.5,0,0), oma=c(0.5,0.5,0.5,0.5), mgp = c(1.5,0.5,0))
plot(NMDS, type = "n")
points(NMDS, pch = c(16,17)[localf], col = c('black', 'grey50')[localf])
ordihull(NMDS, envs$local, lty = c('solid', 'dotted'))
plot(fit_env, col = "grey20")
mtext(paste("Stress:", round(NMDS$stress, 3), " "), side = 3, adj = 1, line = -1.2)
legend("topleft", legend=c('North','South'), col=c('black','grey50'), lty=c('solid', 'dotted'), pch = c(16,17), cex=1)
dev.off()




# PLOT - riqueza x local
# png_fit("rarefacao2.png", pt = 8, width = "14cm", height = "8.6cm", res= 1000)
# par(mar=c(2.5,2.5,0,0), oma=c(0.5,0.5,0.5,0.5), mgp = c(1.5,0.5,0))
set.seed(12345)
xsul = specaccum(mc[envs$local == "South",], method = "random")
plot(xsul, col = 'grey50', ci.type = "polygon", ci.lty = 0, ci.col = 'grey90', cex.lab=0.8, cex.axis=0.7,
     ylim=c(0,120), xlim=c(1,20), xaxs = "i", yaxs = "i", xlab = "Randomized Dredges", ylab="Richness")

xnorte = specaccum(mc[envs$local == "North",], method = "random")
plot(xnorte, col = 'grey20', ci.type = "polygon", ci.lty = 0, ci.col = adjustcolor('grey30', alpha.f = 0.5), add = T)
legend("topleft", legend=c('North','South'), col=c('grey30','grey70'), lty=1, cex=1)
# dev.off()







# output Table 1
envs = dados[, c(.(Ninds = sum(Ninds)), lapply(.SD, first)), by = .(cruzeiro, operacao), .SDcols = c('lon', 'lat', 'shape_leng', 'depth_gebco', 'depth', 'crusts')]
setorder(envs, cruzeiro, operacao)

mc = dcast(dados, cruzeiro + operacao ~ morfotipo, value.var = 'Ninds', fun = sum)
setorder(mc, cruzeiro, operacao)
mc[, c('cruzeiro', 'operacao') := NULL]

envs$S <- specnumber(mc)
envs$H <- round(diversity(mc, index = "shannon"), 3)
envs$J <- round(envs$H / log(envs$S), 3)



envs[is.na(depth), depth := depth_gebco]
envs[, depth := round(depth) * -1]
envs[, depth_gebco := NULL]

envs[, shape_leng := round(shape_leng)]

gd2gm <- function(coord, suffix) {
  coord = abs(coord)
  g = floor(coord)
  sprintf("%d°%04.1f'%s", g, (coord - g)*60, suffix)
}

envs[, lon := gd2gm(lon, "W")]
envs[, lat := gd2gm(lat, "S")]

setcolorder(envs, c(1,2,5,4,6,7,8,3,9:11))
colnames(envs) = c('Cruise', 'Dredge Number', 'Latitude', 'Longitude', 'Length (m)', 'Depth (m)', 'Fe-Mn Crust coverage (%)', 'No. inds.', 'S', "H'", "J'")

fwrite(envs, "main_table.csv")



# PLOT - north vs south relative abundance
tmp = dados2[, .(Ninds = sum(Ninds)), by = .(local, grupo)]
g <- ggplot(tmp, aes(x=Ninds, y=local)) +
  geom_bar(aes(fill=grupo), stat = 'identity', position=position_fill(reverse = TRUE), color = 'black', size = 0.2) + 
  scale_fill_brewer(palette = "Greys", direction = 1) +
  scale_x_continuous(expand = c(0,0), labels = scales::percent) +
  labs(fill = "Phylum", pattern = "Phylum", x = "Relative abundance (%)") +
  theme_classic(base_size = 9, base_line_size = 0.25) +
  theme(axis.title.y = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(size = rel(0.9)),
        legend.key.size = unit(3.5, "mm"))

g

ggsave("north_south.png", g, height = 4, width = 12, units = 'cm', dpi = 600)






# output Table 1

hosp[grupo == "outros", grupo := "Others"]

temp = hosp[, .N, by = .(hosp_morfotipo, hosp_taxon)]
temp[, hosp_taxon := make.names(hosp_taxon, unique = TRUE)]
hosp[temp, hosp_taxon := i.hosp_taxon, on = .(hosp_morfotipo)]



concat <- function(x) {
  xx = table(x)
  paste(names(xx), xx, collapse = ", ")
}

out <- hosp[, .(Ninds = sum(Ninds), S = length(unique(morfotipo)), list = concat(hosp_taxon)), by = .(grupo)]
fwrite(out, "hosp.csv")




# GLM models

fit = glm(Ninds ~ local + dist_rift + crusts + depth + slope + mean_dir, quasipoisson, envs)
summary(fit)
plot(fit)
anova(fit, test= "F")


fit = glm(S ~ local + dist_rift + crusts + depth + slope + mean_dir, quasipoisson, envs)
summary(fit)
plot(fit)


envs[, .(Nidns = mean(Ninds)), by=.(local)]

plot(Ninds ~ crusts, envs)
