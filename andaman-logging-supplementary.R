
rm(list = ls())
library(raster)
library(sf)
library(tidyverse)

#library(rgbif)
#library(Rfast)
# plotcenter <-
#   read_csv("./inputfiles/plotcentroids.csv")
# db_geo <-
#   plotcenter %>%
#   dplyr::mutate(decimalLatitude = latitude,decimalLongitude = longitude) %>%
#   dplyr::select(-latitude,-longitude,-treat) %>%
#   elevation(username = "akshayS",elevation_model = "srtm1") #30x30, register on elevation rgbif API with your username to activate
#
# srtm_ani <-
#   raster("D:/Resources_Grants_Protocols_Datasets_Stats_Teach/data/srtm_ani/N12E092.hgt")
#
# srtm_aniUTM <- projectRaster(from = srtm_ani,
#                              crs = crstext_utm46N)
# slope_raster <- terrain(x = srtm_aniUTM,
#                         opt = "slope",
#                         unit = "degrees",
#                         neighbors = 4)
# db_geo$slope <- extract(x = slope_raster, y = plotcenter_UTM)
#
# plot_pairdist <- st_distance(plotcenter_UTM)
# colnames(plot_pairdist) <- sort(db_geo$plot_ID)
#
# db_fgeo <-
# left_join(x = db_geo,
#           y = data.frame(
#             plotDis_nn_m = apply(X = plot_pairdist,
#                              MARGIN = 2,
#                              FUN = function(X) Rfast::nth(x = X,k = 2))) %>%
#             rownames_to_column("plot_ID"),
#           by = "plot_ID") %>%
#   rename(elev = elevation_geonames) %>%
#   mutate(elev = ifelse(elev<0,0,elev))
# write_csv(x = db_fgeo,path = "db_fgeo.csv")

# db <-
#   read_csv("./inputfiles/plotdata.csv") %>%
#   mutate(treat = case_when(treatment ==      "TL" ~ "L2",
#                            treatment == "OL_near" ~ "L1",
#                            treatment ==  "OL_far" ~ "B",
#                            treatment =="baseline" ~ "B"),
#          treat = factor(treat, levels = c("B","L1","L2")),
#          forest_type = factor(forest_type, levels = c("deciduous","evergreen")),
#          plot_ID = factor(plot_ID)) %>%
#   dplyr::select(-treatment)
#
# db_fgeo <- read_csv("db_fgeo.csv") %>%
#   left_join(x = .,y = db %>% dplyr::select(plot_ID,treat) %>% distinct())
#
# plotcenter_UTM  <-
#   read_sf("./inputfiles/plotcentroidsUTM64N.shp")
#
# crstext_utm46N <- crs(plotcenter_UTM)
# #'+proj=utm +zone=46 +datum=WGS84 +units=m +no_defs '
#
# resettle <- matrix(data = NA,ncol = 2,nrow = 11)
# colnames(resettle) <- c("is","becomes")
# resettle[,"is"]      <- c(1,2,3, 5, 7,8,9,10,11,13,16)
# resettle[,"becomes"] <- c(1,1,1,NA,NA,NA,NA,NA,NA,NA,NA)
#
# dd <- raster(x = "inputfiles/LULC_Bhuvan_studyregion.tif")
# dd_settlement <- reclassify(x = dd,rcl = resettle)

# dd_settlement2 <-
#   rasterToPolygons(x = dd_settlement,n = 4,na.rm = T,digits = 5)
# dd_settlement3 <- sf::st_as_sf(dd_settlement2)
# dd_settlement4 <- st_transform(x = dd_settlement3,
#                                crs = crstext_utm46N)
#
# distmat <- st_distance(x = plotcenter_UTM,y = dd_settlement4)
#
# distance_to_settlement <- apply(X = distmat,MARGIN = 1,FUN = min)
# db_fgeo$dist_settlement <- distance_to_settlement
#
# write_csv(x = db_fgeo,path = "October 2020 - Major revision/FigureS0_inputtable.csv")
#
# aniext <- as(raster::extent(92.70, 92.90, 12.10, 12.40), "SpatialPolygons")
# #extent of 50 plots #in SL
# proj4string(aniext) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
#
# temp0 <- raster(x = "E:/WorldClimData_200GB/wc2.1_30s_bio/wc2.1_30s_bio_1.tif")
# # annual mean temperature
# rain0 <- raster(x = "E:/WorldClimData_200GB/wc2.1_30s_bio/wc2.1_30s_bio_12.tif")
# # annual precipitation
#
# temp <- crop(x = temp0,y = aniext)
# rain <- crop(x = rain0,y = aniext)
#
# coords <-
#   db_fgeo[,c("latitude","longitude")] %>%
#   st_as_sf(x = .,coords = c("longitude","latitude"),
#            crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
#
# db_fgeo$temperature <- raster::extract(temp,coords)
# db_fgeo$rainfall <- raster::extract(rain,coords)
#
# write_csv(x = db_fgeo,path = "October 2020 - Major revision/FigureS0_inputtable.csv")
# write_csv(x = db_fgeo,path = "db_fgeo.csv")

# v1 <-
# db_fgeo %>%
#   group_by(treat,forest_typ) %>%
#   summarise(median_val = median(slope),
#             min_val = min(slope), max_val = max(slope))

# repeat above for each variable: elevation, slope, distance_edge, distance_nearestneighbour


#### Fig S2 - Identifying and removing plots that are erroneously classified ####
# (see supplementary.R script)

## k-means clustering (k=2) on compositional difference
## n=4 plots with wrong ID
## NDMS of all plots (P.S. k=3 in metaMDS() corresponds to 3 dof, not 3 clusters).
## coloured according to user-assigned forest type, to visually confirm misclassification

detach(name = "package:raster",unload = TRUE)
db0 <-
     read_csv("plotdata.csv") %>%
     mutate(treat = case_when(treatment ==      "TL" ~ "L2",
                              treatment == "OL_near" ~ "L1",
                              treatment ==  "OL_far" ~ "B",
                              treatment =="baseline" ~ "B"),
            treat = factor(treat, levels = c("B","L1","L2")),
            forest_type = factor(forest_type, levels = c("deciduous","evergreen")),
            plot_ID = factor(plot_ID)) %>%
     dplyr::select(-treatment)

nmdsdat1 <-
  db0 %>%
  filter(gbh_cm>=30,
         plotA=="Y",
         !code %in% c("stump","belsp","unk"),
         gf!="P",
         ID_level!="F") %>%
  select(forest_type,plot_ID,code,plotA, plotC)
# 69% of all rows

dbsup <-
  dcast(data = nmdsdat1, formula = plot_ID + forest_type ~ code,
        fun.aggregate = length, value.var = "code") %>%
  column_to_rownames("plot_ID")

dbsup2 <- decostand(x = dbsup %>% select(-forest_type), method = "hellinger")

nmdsdat2 <- metaMDS(comm =  dbsup2,
                    distance = "bray",
                    k = 2,try = 300,trymax = 1000)

dbsup2.nmds <-
  nmdsdat2$points[,c(1:2)] %>%
  cbind() %>%
  as.data.frame()

if(identical(row.names(dbsup2.nmds),row.names(dbsup2)))
{   dbsup2.nmds$forest_type<-dbsup2$forest_type } #just a check

tmp3 <- kmeans(x = dbsup2, centers = 2)
tmp4 <-
  left_join(x = tmp3$cluster %>% as.matrix() %>% as.data.frame() %>% rownames_to_column(),
            y = nmdsdat1 %>% group_by(plot_ID, forest_type) %>% summarise(nn = n()) %>% select(-nn),
            by = c("rowname" = "plot_ID")) %>%
  mutate(tag = paste0(V1,"_",forest_type))

E.code <- ifelse(test = length(which(tmp4$tag=="1_evergreen")) > length(which(tmp4$tag=="2_evergreen")),
                 yes = 1,no = 2)

D.code <-ifelse(test = E.code==1,yes = 2,no = 1)

removeplots <- tmp4[with(tmp4,which((V1==E.code & forest_type!="evergreen") | (V1==D.code & forest_type!="deciduous"))),'rowname']

dbsup2.nmds <-
  dbsup2.nmds %>%
  rownames_to_column("plot_ID")

beforeremplot <- ggplot() +
  geom_point(data = left_join(x = dbsup2.nmds,y = db0 %>% select(forest_type,plot_ID) %>% distinct()),
             aes(y = MDS2,
                 x = MDS1,
                 colour = forest_type),
             size = 3) +
  geom_point(data = dbsup2.nmds %>% filter(plot_ID %in% removeplots),
             aes(y=MDS2,x=MDS1), colour="black", size = 4) +
  theme_light() + ggtitle("") +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5)) +
  scale_colour_manual(name = "", values = c("dodgerblue2","forestgreen")) + xlim(-1.25,1) + ylim(-1.25,1)

ggsave(plot = beforeremplot,
       filename = "output/figS1_NMDS.png",
       device = "png", width = 5.5, height = 5.5,units = "in", dpi = 300)


#### Table S1 - Dominant species ####
tmp1 <-
  db0 %>%
  select(gbh_cm,gbh_cm_stem2,gbh_cm_stem3) %>%
  replace_na(list(gbh_cm_stem2=0,
                  gbh_cm_stem3=0))

db0$girth_eff <-
  apply(X = tmp1,MARGIN = 2,FUN = function(x) x^2) %>%
  rowSums() %>% sqrt()

cutoff<-0.5 #50% of all stems
samplingscaling<-c(12,12,18,11,12,7)
dbdom<-db0 %>% filter(plotA=="Y" | plotB=="Y") #only in treefall plot

############# adult tree dominants
dom_medtree<-
  db0 %>%
  filter(girth_eff>=30,code!="stump",code!="belsp",gf!="P") %>%
  group_by(treat,forest_type,species_ID) %>%
  summarise(ntree = n()) %>%
  dcast(species_ID ~ treat + forest_type, value.var = "ntree") %>%
  replace_na(replace = list(B_deciduous=0,B_evergreen=0,
                            L1_deciduous=0,L1_evergreen=0,
                            L2_deciduous=0,L2_evergreen=0)) %>%
  column_to_rownames("species_ID")

dom_medtree1 <-
  round(t(t(dom_medtree)/samplingscaling),2) %>%
  as.data.frame() %>% rownames_to_column("species_ID")

b_d <- dom_medtree1 %>% select(species_ID,B_deciduous) %>% arrange(-B_deciduous) %>% mutate(cums = cumsum(B_deciduous),props = cums/max(cums)) %>% filter(props<=cutoff) %>% select(species_ID)
l1_d <- dom_medtree1 %>% select(species_ID,L1_deciduous) %>% arrange(-L1_deciduous) %>% mutate(cums = cumsum(L1_deciduous),props = cums/max(cums)) %>% filter(props<=cutoff) %>% select(species_ID)
l2_d <- dom_medtree1 %>% select(species_ID,L2_deciduous) %>% arrange(-L2_deciduous) %>% mutate(cums = cumsum(L2_deciduous),props = cums/max(cums)) %>% filter(props<=cutoff) %>% select(species_ID)
b_e <- dom_medtree1 %>% select(species_ID,B_evergreen) %>% arrange(-B_evergreen) %>% mutate(cums = cumsum(B_evergreen),props = cums/max(cums)) %>% filter(props<=cutoff) %>% select(species_ID)
l1_e <- dom_medtree1 %>% select(species_ID,L1_evergreen) %>% arrange(-L1_evergreen) %>% mutate(cums = cumsum(L1_evergreen),props = cums/max(cums)) %>% filter(props<=cutoff) %>% select(species_ID)
l2_e <- dom_medtree1 %>% select(species_ID,L2_evergreen) %>% arrange(-L2_evergreen) %>% mutate(cums = cumsum(L2_evergreen),props = cums/max(cums)) %>% filter(props<=cutoff) %>% select(species_ID)
med <- Reduce(union,list(b_d,l1_d,l2_d,b_e,l1_e,l2_e)) #25 species

table1_adult<-
  dom_medtree1 %>% filter(species_ID %in% med$species_ID) %>% select(species_ID,L2_deciduous,L1_deciduous,B_deciduous,L2_evergreen,L1_evergreen,B_evergreen)

#############  pole tree dominants #
dom_smatree<-
  db0 %>% filter(girth_eff<30) %>% filter(code!="stump",code!="belsp",gf!="P") %>% group_by(treat,forest_type,species_ID) %>% summarise(ntree = n()) %>% dcast(species_ID ~ treat + forest_type, value.var = "ntree") %>% replace_na(replace = list(B_deciduous=0,B_evergreen=0,L1_deciduous=0,L1_evergreen=0,L2_deciduous=0,L2_evergreen=0)) %>% column_to_rownames("species_ID")
samplingscaling<-c(12,12,18,11,12,7)
dom_smatree1<-round(t(t(dom_smatree*5)/samplingscaling),2) %>% as.data.frame() %>% rownames_to_column("species_ID") #scaled to 0.05 ha
b_d<-dom_smatree1 %>% select(species_ID,B_deciduous) %>% arrange(-B_deciduous) %>% mutate(cums = cumsum(B_deciduous),props = cums/max(cums)) %>% filter(props<=cutoff) %>% select(species_ID)
l1_d<-dom_smatree1 %>% select(species_ID,L1_deciduous) %>% arrange(-L1_deciduous) %>% mutate(cums = cumsum(L1_deciduous),props = cums/max(cums)) %>% filter(props<=cutoff) %>% select(species_ID)
l2_d<-dom_smatree1 %>% select(species_ID,L2_deciduous) %>% arrange(-L2_deciduous) %>% mutate(cums = cumsum(L2_deciduous),props = cums/max(cums)) %>% filter(props<=cutoff) %>% select(species_ID)
b_e<-dom_smatree1 %>% select(species_ID,B_evergreen) %>% arrange(-B_evergreen) %>% mutate(cums = cumsum(B_evergreen),props = cums/max(cums)) %>% filter(props<=cutoff) %>% select(species_ID)
l1_e<-dom_smatree1 %>% select(species_ID,L1_evergreen) %>% arrange(-L1_evergreen) %>% mutate(cums = cumsum(L1_evergreen),props = cums/max(cums)) %>% filter(props<=cutoff) %>% select(species_ID)
l2_e<-dom_smatree1 %>% select(species_ID,L2_evergreen) %>% arrange(-L2_evergreen) %>% mutate(cums = cumsum(L2_evergreen),props = cums/max(cums)) %>% filter(props<=cutoff) %>% select(species_ID)
sma<-Reduce(union,list(b_d,l1_d,l2_d,b_e,l1_e,l2_e)) #125 species

table1_pole<-
  dom_smatree1 %>%
  filter(species_ID %in% sma$species_ID) %>%
  select(species_ID,L2_deciduous,L1_deciduous,B_deciduous,L2_evergreen,L1_evergreen,B_evergreen)

names(table1_pole) <- c("speciesID","L2(D,p)","L1(D,p)","B(D,p)","L2(E,p)","L1(E,p)","B(E,p)")
names(table1_adult) <- c("speciesID","L2(D,a)","L1(D,a)","B(D,a)","L2(E,a)","L1(E,a)","B(E,a)")

table1_both <- inner_join(x = table1_pole,table1_adult)
table1_adultsonly <- table1_adult %>% filter(!speciesID %in% table1_both$speciesID)
table1_poleonly <- table1_pole %>% filter(!speciesID %in% table1_both$speciesID)

sink("output/tableS1_andaman-logging.csv")
print("Species that appeared as both pole-sized and adult-sized trees only\n")
table1_both
print("Species that appeared as pole-sized trees only\n")
table1_poleonly
print("Species that appeared as adult-sized trees only\n")
table1_adultsonly
sink()









