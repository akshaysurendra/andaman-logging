## making a map (figure 1 - A,B,C)
library(tidyverse)
library(sf)
library(cowplot)
library(OpenImageR)
plotloc <- read_csv("inputfiles/plotcentroids.csv")
indiast <- read_sf("inputfiles/IndiaStates_SHP/INDIA.shp")
ani <- indiast %>% filter(ST_NAME=="ANDAMAN AND NICOBAR ISLANDS")
ani_bb <- st_bbox(ani)
and_bb <- ani_bb
and_bb[2] <- 10.479123 #ymin brought up to exclude the nicobar islands
and <- 
  st_crop(x = ani,
          y = 
            as(raster::extent(92.20725, 93.5, 10.47913, 13.9000),"SpatialPolygons") %>% 
            st_as_sf() %>% 
            st_set_crs(.,value = st_crs(ani))
  ) %>% 
  mutate(ST_NAME="Andaman Islands") #figure1-B
india <- read_sf("inputfiles/India_boundary_Igismap/India_Boundary.shp") #figure1-A

fig1 <- vector(mode = "list",length = 4)
names(fig1) <- c("A","B","C","D")
fig1$A <- 
  ggplot(india) + geom_sf(fill = "white",colour = "gray30") + 
  geom_rect(xmin = 92,xmax = 94.5, ymin = 10, ymax = 14,
            linejoin = "bevel",linetype = 1, alpha = 0, colour = "gray10", size = 0.8) +
  theme_classic() + theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank()) # + panel_border(colour = "gray50")
ggsave(filename = "inputfiles/fig1a.png",plot = fig1$A,device = "png",width = 1.75,height = 4, dpi = 600, units = "in")
fig1$A <- ggdraw() + draw_image(readImage(path = "inputfiles/fig1a.png"))

fig1$B <- 
  ggplot(and) + geom_sf(fill = "white",colour = "gray30") + 
  geom_rect(xmin = 92.6,xmax = 93.06, ymin = 12, ymax = 12.46,linejoin = "bevel",linetype = 1, alpha = 0, colour = "gray10", size = 0.8) +
  theme_classic() + theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank())# + panel_border(colour = "gray50")
ggsave(filename = "inputfiles/fig1b.png",plot = fig1$B,device = "png",dpi = 600, width = 1.6, height = 4, units = "in")
fig1$B <- ggdraw() + draw_image(readImage(path = "inputfiles/fig1b.png"))

tmp_c <- readImage(path = "inputfiles/fig1C_fin_unscaled.png")
ggsave(filename = "inputfiles/fig1C_fin.png",plot = ggdraw() + draw_image(tmp_c),device = "png",dpi = 600, width = 4.8, height = 4, units = "in")
fig1$C <- ggdraw() + draw_image(readImage(path = "inputfiles/fig1C_fin.png"))

tmp_d <- readImage(path = "inputfiles/fig1D.png")
ggsave(filename = "inputfiles/fig1D_fin.png",plot = ggdraw() + draw_image(tmp_d),device = "png",dpi = 600, width = 8, height = 4, units = "in")
fig1$D <- ggdraw() + draw_image(readImage(path = "inputfiles/fig1D_fin.png"))


toprow <- plot_grid(plotlist = fig1[c(1,2,3)],nrow = 1,ncol = 3,rel_widths = c(1,1,3),labels = "AUTO",label_size = 16,greedy = TRUE, align = "h")
ggsave(filename = "inputfiles/toprow.png",plot = toprow,width = 8, height = 4,units = "in", dpi = 600)

bottomrow <- plot_grid(plotlist = fig1[c(4)],nrow = 1, ncol = 1, labels = "D",label_size = 16,greedy = TRUE, align = "h")
ggsave(filename = "inputfiles/bottomrow.png",plot = bottomrow,width = 8, height = 4,units = "in", dpi = 600)
fig1fin <- plot_grid(plotlist = list(toprow,bottomrow),nrow = 2,ncol = 1,greedy = TRUE)
ggsave(filename = "inputfiles/fig1_FIN.png",plot = fig1fin,scale = 2,width = 8, height = 8, units = "in", dpi = 600)


#write_sf(and,"andamanoutline.shp")






dd <- raster(x = "inputfiles/LULC_Bhuvan_studyregion.tif")

reclass <- matrix(data = NA,ncol = 2,nrow = 11)
colnames(reclass) <- c("is","becomes")
reclass[,"is"]      <- c(1,2,3, 5, 7,8,9,10,11,13,16)
reclass[,"becomes"] <- c(1,1,1,NA,NA,2,2, 2, 3, 3,NA)

dd_new <- reclassify(x = dd,rcl = reclass)



# writeRaster(x = dd_new,filename = "newLULC_reclass.tif",overwrite=TRUE)
#ggplot(data = dd_new) + geom_sf() + annotate_scal

# studysite <- rasterToPolygons(x = dd_new,na.rm = TRUE,dissolve = TRUE)
#studysite <- polygonize(x = dd_new)
#write_sf(studysite,"studysite_polygon.shp")



p1 <- 
  gplot(dd_new) + geom_tile(aes(fill=factor(value, labels = c("Habitation","Forest","Mangrove")))) + 
  scale_fill_manual(values = c("grey10","grey50","grey80"), name = "") + theme_classic() +
  geom_polygon(data = and,aes(x=Long,y=Lat)) +
theme(legend.position = "top") + xlab("") + ylab("")
ggsave(filename = "maptrial_1C.png",plot = p1,device = "png",width = 14,height = 14,units = "cm",dpi = 600)


