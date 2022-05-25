library(sf)
library(raster)
library(data.table)

# Get variables from rasters (depth, slope, mean direction)
# Used for generate_table.R

# check N.B. Morgan et al. / Deep-Sea Research I 104 (2015) 92-105
mean_dir <- function(C, S) {
  x1=atan(S/C)
  x2=ifelse(S > 0 && C > 0, 0, ifelse(C < 0, pi, 2*pi))
  return(x1 + x2)
}


shp <- list(rgr1 = read_sf(dsn = "G:\\Meu Drive\\Doutorado\\cruzeiro RGR1", layer = "dredge"),
            dy94 = read_sf(dsn = "G:\\Meu Drive\\Doutorado\\cruzeiro DY94", layer = "dredges"),
            hybis = read_sf(dsn = "G:\\Meu Drive\\Doutorado\\cruzeiro DY94", layer = "hybis")
)

op <- c('N_operaç', 'operacao', 'dive')
cruzeiro <- c('RGR1', 'DY94', 'DY94')

# load rasters
rasters = c("depth", "slope", "cos_aspect", "sin_aspect")
rasters = paste0("G:\\Meu Drive\\Doutorado\\images and data\\BTM_DY94\\", rasters, ".tif")
r = stack(rasters)

gebco = raster("G:\\Meu Drive\\ArcGis\\gebco_08a_Subset0.img")
            
out <- list()
for (i in 1:length(shp)) {
  temp = st_transform(shp[[i]], CRS("+proj=utm +zone=25 +south +datum=WGS84 +units=m +no_defs"))
  
  
  # get points every X meters
  temp2 = st_line_sample(temp, density = 1/5)
  temp2 = st_coordinates(temp2)
  
  x = extract(r, temp2[, 1:2])
  x = data.table(operacao = temp[[op[i]]][temp2[, 3]], x)
  
  # get mean of variables
  x1 = x[, lapply(.SD, mean), by = .(operacao), .SDcols = c("depth", "slope")]
  x2 = x[, lapply(.SD, sum), by = .(operacao), .SDcols = c("cos_aspect", "sin_aspect")]
  
  x1$mean_dir = mean_dir(x2$cos_aspect, x2$sin_aspect) * 180 / pi
  x1$cruzeiro = cruzeiro[i]
  x1$lon = temp$MID_X
  x1$lat = temp$MID_Y
  x1$depth_gebco = extract(gebco, x1[, .(lon, lat)])
  
  out[[i]] = copy(x1)
}


out = rbindlist(out)
out[, operacao := sub("D([0-9])$", "D0\\1", operacao)]

fwrite(out, "draga_vars.csv")
