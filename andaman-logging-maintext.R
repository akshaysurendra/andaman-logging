# Impact of repeat logging on forest recovery
# Akshay Surendra, Anand M Osuri, Jayashree Ratnam
# August 2020

#### Packages, data input, data prep ######
#setwd("~/Desktop/draft/FEM manuscript/")

rm(list = ls())
library(vegan)
library(tidyverse)
library(gridExtra)
library(reshape2)
library(grid)
library(pander)
library(ggsignif)
library(ggmap)
library(cowplot)
library(betareg)
library(sjPlot)


### Table 1 #####

db_fgeo <- read_csv("db_fgeo.csv")

tab1 <- 
  full_join(x = db_fgeo %>% group_by(treat,forest_typ) %>% summarise(Nplots = n_distinct(plot_ID)),
            y = db_fgeo %>% 
              group_by(treat,forest_typ) %>% 
              summarise_at(.vars = c("temperature","rainfall","elev","slope","dist_settlement"),
                           .funs = function(x) paste0(round(median(x,na.rm=T),1),
                                                      " (",round(min(x,na.rm=T),1),
                                                      "-",round(max(x,na.rm=T),1),
                                                      ")")
                           ),
            by = c("forest_typ","treat"))

# write_csv(x = tab1,path = "October 2020 - Major revision/Table1_plotchar.csv")


#### ######
db0 <- read.csv("./inputfiles/plotdata.csv",header = T, stringsAsFactors = F)

### Identifying and removing plots that are erroneously classified
# (see supplementary.R script)
removeplots<-c("BANL0303","BAT1P10","MAT1C6P09","MAT2C2P09")

dbtmp0 <- 
  db0 %>% 
  filter(!plot_ID %in% removeplots) #### 4 plots removed


### Re-level treatment, compute AGB, join phenology (D/E) & canopy data ####
  
dbtmp <- 
  dbtmp0 %>% 
  mutate(treat = case_when(treatment ==      "TL" ~ "L2",
                           treatment == "OL_near" ~ "L1",
                           treatment ==  "OL_far" ~ "B",
                           treatment =="baseline" ~ "B"),
         treat = factor(treat, levels = c("B","L1","L2")),
         forest_type = factor(forest_type, levels = c("deciduous","evergreen")),
         plot_ID = factor(plot_ID)) %>% 
  select(-treatment)

dbtmp1 <- 
  dbtmp %>% 
  mutate(agb = data.frame(agbstem1 = ( (0.0673) * ((dbtmp$wd_m) * (dbtmp$dia1^2) * (dbtmp$height_m))^(0.976) ),
                          agbstem2 = ( (0.0673) * ((dbtmp$wd_m) * (dbtmp$dia2^2) * (dbtmp$height_m))^(0.976) ),
                          agbstem3 = ( (0.0673) * ((dbtmp$wd_m) * (dbtmp$dia3^2) * (dbtmp$height_m))^(0.976) )) %>% 
           apply(X = ., MARGIN = 1, FUN = sum, na.rm = T), 
         agc_Mg = 0.5*agb/1000)
#Chave 2014 pan-tropical formula used; Carbon fraction ~50% of AGB (Chave 2005)

phenology <- read_csv("inputfiles/treecodes.csv")
db <-
  dbtmp1 %>% 
  left_join(x = .,y = phenology %>% select(species_ID,phenologyFIN))
setdiff(x = unique(db$species_ID),
        y = unique(phenology$species_ID)) # only stump has no phenology


cc <- read.csv("./inputfiles/canopy.csv",header = T, stringsAsFactors = F)
cc$lowerstorey <- apply(X = cc[,c(2:10)],
                         MARGIN = 1,FUN = sum)
cc$upperstorey <- apply(X = cc[,c(11:19)],
                         MARGIN = 1,FUN = sum)
## total canopy over 45 metres summed up by adding counts from 9 5-metre sections

canopy <- 
  cc %>% 
  rename(plot_ID = "PLOT") %>% 
  select(plot_ID,lowerstorey, upperstorey) %>% 
  mutate(plot_ID = factor(plot_ID),
         lowerstorey = 100*lowerstorey/45,
         upperstorey = 100*upperstorey/45) %>% 
  right_join(x = .,
             y = db %>% select(plot_ID,forest_type,treat) %>% distinct()) %>% 
  # right join as 4 plots already removed from db but not cc
  pivot_longer(cols = c("lowerstorey","upperstorey"),
               names_to = "strata",values_to = "ccover") %>% 
  mutate(ccover = ccover/100)
  # percentage to proportion

cbPalette <- rev(c("#999999","#56B4E9", "#F0E442","#D55E00"))

dbp <- modop <- modelF <- modeldf <- vector(mode = "list",length = 18)
dbpsumm <-   vector(mode = "list",length = 18)
names(dbp) <- names(dbpsumm) <- names(modop) <- names(modelF) <- names(modeldf) <- 
  c("stem_pole_deciduous","stem_pole_evergreen",
    "stem_adult_deciduous","stem_adult_evergreen",
    "understory_deciduous","understory_evergreen",
    "overstory_deciduous","overstory_evergreen",
    "richness_pole_deciduous","richness_pole_evergreen",
    "richness_adult_deciduous","richness_adult_evergreen",
    "decifrac_pole_deciduous","decifrac_pole_evergreen",
    "decifrac_adult_deciduous","decifrac_adult_evergreen",
    "AGB_deciduous","AGB_evergreen")

      dbstem <- db %>% filter(code != "stump", code != "belsp", gf!="P") # identity level doesn't matter, no palms, lianas, stumps
 dbstem.pole <- dbstem %>% filter(gbh_cm>=10, gbh_cm <30) #pole-sized trees
dbstem.adult <- dbstem %>% filter(gbh_cm>=30, plotA=="Y") #adults

dbspec <- db %>% 
  filter(!code %in% c("stump","belp","unk"), gf!="P", ID_level!="F")
dbspec.pole <- dbspec %>% filter(gbh_cm>=10, gbh_cm<30)
dbspec.adult <- dbspec %>% filter(gbh_cm>=30, plotA=="Y")

dbagc.adult <-
  db %>%
  filter(code!="stump", code!="belsp", gf!="P", gbh_cm>=30, plotA=="Y")


#### Stem density ~ logging frequency (pole) ####

dbp$stem_pole_deciduous <- 
  dbstem.pole %>% 
  filter(forest_type=="deciduous") %>% 
  group_by(treat,plot_ID) %>% summarise(nstem = n())
dbpsumm$stem_pole_deciduous <- 
  dbp$stem_pole_deciduous %>% group_by(treat) %>%
  summarise(y = mean(nstem), ns = n(), ysd = sd(nstem)) %>% 
  mutate(ymin = y - 1.96*(ysd/sqrt(ns)), ymax = y + 1.96*(ysd/sqrt(ns)),
         forest_type = "deciduous")

dbp$stem_pole_evergreen <- 
  dbstem.pole %>% 
  filter(forest_type=="evergreen") %>% 
  group_by(treat,plot_ID) %>% summarise(nstem = n()) %>% 
  ungroup() %>% 
  add_row(treat = "B", plot_ID = "BAT2P13", nstem = 0) #0 stems so missed out
dbpsumm$stem_pole_evergreen <- 
  dbp$stem_pole_evergreen %>% group_by(treat) %>%
  summarise(y = mean(nstem), ns = n(), ysd = sd(nstem)) %>% 
  mutate(ymin = y - 1.96*(ysd/sqrt(ns)), ymax = y + 1.96*(ysd/sqrt(ns)),
         forest_type = "evergreen")

#### stem density ~ logging frequency (adult) ####

dbp$stem_adult_deciduous<-dbstem.adult %>% 
  filter(forest_type=="deciduous") %>%
  group_by(treat,plot_ID) %>% summarise(nstem = n())
dbpsumm$stem_adult_deciduous <- 
  dbp$stem_adult_deciduous %>% group_by(treat) %>%
  summarise(y = mean(nstem), ns = n(), ysd = sd(nstem)) %>% 
  mutate(ymin = y - 1.96*(ysd/sqrt(ns)), ymax = y + 1.96*(ysd/sqrt(ns)),
         forest_type = "deciduous")

dbp$stem_adult_evergreen<-dbstem.adult %>%
  filter(forest_type=="evergreen") %>%
  group_by(treat,plot_ID) %>% summarise(nstem = n())
dbpsumm$stem_adult_evergreen <- 
  dbp$stem_adult_evergreen %>% group_by(treat) %>%
  summarise(y = mean(nstem), ns = n(), ysd = sd(nstem)) %>% 
  mutate(ymin = y - 1.96*(ysd/sqrt(ns)), ymax = y + 1.96*(ysd/sqrt(ns)),
         forest_type = "evergreen")

#### canopy cover ~ logging frequency (understory) ####

dbp$understory_deciduous <- canopy %>% filter(strata =="lowerstorey",
                                              forest_type=="deciduous")
dbpsumm$understory_deciduous <- 
  dbp$understory_deciduous %>%
  filter(forest_type=="deciduous",strata=="lowerstorey") %>% 
  group_by(treat) %>% 
  summarise(ns = n(),
            y = mean(ccover), sd = sd(ccover),
            ymax = y + 1.96*sd/sqrt(ns), ymin = y - 1.96*sd/sqrt(ns),
            forest_type = "deciduous")

dbp$understory_evergreen <- canopy %>% filter(strata =="lowerstorey",
                                              forest_type=="evergreen")
dbpsumm$understory_evergreen <- 
  dbp$understory_evergreen %>%
  filter(forest_type=="evergreen",strata=="lowerstorey") %>% 
  group_by(treat) %>% 
  summarise(ns = n(),
            y = mean(ccover), sd = sd(ccover),
            ymax = y + 1.96*sd/sqrt(ns), ymin = y - 1.96*sd/sqrt(ns),
            forest_type = "evergreen")

dbp$overstory_deciduous <- canopy %>% filter(strata =="upperstorey",
                                             forest_type=="deciduous")
dbpsumm$overstory_deciduous <- 
  dbp$overstory_deciduous %>%
  filter(forest_type=="deciduous",strata=="upperstorey") %>% 
  group_by(treat) %>% 
  summarise(ns = n(),
            y = mean(ccover), sd = sd(ccover),
            ymax = y + 1.96*sd/sqrt(ns), ymin = y - 1.96*sd/sqrt(ns),
            forest_type = "deciduous")

dbp$overstory_evergreen <- canopy %>% filter(strata =="upperstorey",
                                             forest_type=="evergreen")
dbpsumm$overstory_evergreen <- 
  dbp$overstory_evergreen %>%
  filter(forest_type=="evergreen",strata=="upperstorey") %>% 
  group_by(treat) %>%
  summarise(ns = n(),
            y = mean(ccover), sd = sd(ccover),
            ymax = y + 1.96*sd/sqrt(ns), ymin = y - 1.96*sd/sqrt(ns),
            forest_type = "evergreen")

#### species richness (exp-Shannon) ~ logging frequency (pole) ####

tmpP <- 
  dbspec.pole %>% 
  dcast(data = .,
        formula = plot_ID ~ code,fun.aggregate = length, value.var = "code") %>%
  column_to_rownames("plot_ID") %>%
  diversity(index = "shannon") %>%
  exp() %>%
  data.frame(expH = .) %>% 
  rownames_to_column("plot_ID") %>% 
  add_row(plot_ID = "BAT2P13",expH = 0) %>% 
  left_join(x = .,y = db %>% select(forest_type,treat,plot_ID) %>% distinct())

dbp$richness_pole_deciduous <- tmpP %>% filter(forest_type=="deciduous")
dbpsumm$richness_pole_deciduous <- 
  dbp$richness_pole_deciduous %>%
  group_by(treat) %>% 
  summarise(y = mean(expH), ns = n(), ysd = sd(expH)) %>% 
  mutate(ymin = y - 1.96*(ysd/sqrt(ns)), ymax = y + 1.96*(ysd/sqrt(ns)),
         forest_type = "deciduous")

dbp$richness_pole_evergreen <- tmpP %>% filter(forest_type=="evergreen")
dbpsumm$richness_pole_evergreen<- 
  dbp$richness_pole_evergreen %>%
  group_by(treat) %>%
  summarise(y = mean(expH), ns = n(), ysd = sd(expH)) %>%
  mutate(ymin = y - 1.96*(ysd/sqrt(ns)), ymax = y + 1.96*(ysd/sqrt(ns)),
         forest_type = "evergreen")


#### species richness (exp-Shannon) ~ logging frequency (adult) ####

tmpA <- 
  dbspec.adult %>% 
  dcast(data = .,
        formula = plot_ID ~ code,fun.aggregate = length, value.var = "code") %>%
  column_to_rownames("plot_ID") %>%
  diversity(index = "shannon") %>%
  exp() %>%
  data.frame(expH = .) %>% 
  rownames_to_column("plot_ID") %>% 
  left_join(x = .,y = db %>% select(forest_type,treat,plot_ID) %>% distinct())

dbp$richness_adult_deciduous <- tmpA %>% filter(forest_type=="deciduous")
dbpsumm$richness_adult_deciduous <- 
  dbp$richness_adult_deciduous %>%
  group_by(treat) %>%
  summarise(y = mean(expH), ns = n(), ysd = sd(expH)) %>%
  mutate(ymin = y - 1.96*(ysd/sqrt(ns)), ymax = y + 1.96*(ysd/sqrt(ns)),
         forest_type = "deciduous")

dbp$richness_adult_evergreen <- tmpA %>% filter(forest_type=="evergreen")
dbpsumm$richness_adult_evergreen <-
  dbp$richness_adult_evergreen %>%
  group_by(treat) %>%
  summarise(y = mean(expH), ns = n(), ysd = sd(expH)) %>%
  mutate(ymin = y - 1.96*(ysd/sqrt(ns)), ymax = y + 1.96*(ysd/sqrt(ns)),
         forest_type = "evergreen")

#### deciduous fraction ~  logging frequency (pole) ####

dbp$decifrac_pole_deciduous <- 
  left_join(
    x = dbspec.pole %>% filter(forest_type=="deciduous") %>%
      group_by(treat,plot_ID) %>% summarise(nstem = n()),
    y = dbspec.pole %>% filter(forest_type=="deciduous", phenologyFIN=="D") %>%
      group_by(treat,plot_ID) %>% summarise(decistem = n())) %>% 
  mutate(deciprop = decistem/nstem) %>% 
  replace_na(list(deciprop = 0))

dbpsumm$decifrac_pole_deciduous <- 
  dbp$decifrac_pole_deciduous %>% group_by(treat) %>%
  summarise(y = mean(deciprop), ns = n(), ysd = sd(deciprop)) %>% 
  mutate(ymin = y - 1.96*(ysd/sqrt(ns)), ymax = y + 1.96*(ysd/sqrt(ns)),
         forest_type = "deciduous")

dbp$decifrac_pole_evergreen <- 
  left_join(
    x = dbspec.pole %>% filter(forest_type=="evergreen") %>%
      group_by(treat,plot_ID) %>% summarise(nstem = n()) %>% 
      ungroup() %>% 
      add_row(treat = "B", plot_ID = "BAT2P13", nstem = 0),
    y = dbspec.pole %>% filter(forest_type=="evergreen", phenologyFIN=="D") %>%
      group_by(treat,plot_ID) %>% summarise(decistem = n())) %>% 
  mutate(deciprop = decistem/nstem) %>% 
  replace_na(list(deciprop = 0))

dbpsumm$decifrac_pole_evergreen <- 
  dbp$decifrac_pole_evergreen %>% group_by(treat) %>%
  summarise(y = mean(deciprop), ns = n(), ysd = sd(deciprop)) %>% 
  mutate(ymin = y - 1.96*(ysd/sqrt(ns)), ymax = y + 1.96*(ysd/sqrt(ns)),
         forest_type = "evergreen")

#### deciduous fraction ~  logging frequency (adult) ####

dbp$decifrac_adult_deciduous<- 
  left_join(
    x = dbspec.adult %>% filter(forest_type=="deciduous") %>%
      group_by(treat,plot_ID) %>% summarise(nstem = n()),
    y = dbspec.adult %>% filter(forest_type=="deciduous", phenologyFIN=="D") %>%
      group_by(treat,plot_ID) %>% summarise(decistem = n())) %>% 
  mutate(deciprop = decistem/nstem) %>% 
  replace_na(list(deciprop = 0))

dbpsumm$decifrac_adult_deciduous <- 
  dbp$decifrac_adult_deciduous %>% group_by(treat) %>%
  summarise(y = mean(deciprop), ns = n(), ysd = sd(deciprop)) %>% 
  mutate(ymin = y - 1.96*(ysd/sqrt(ns)), ymax = y + 1.96*(ysd/sqrt(ns)),
         forest_type = "deciduous")

dbp$decifrac_adult_evergreen<- 
  left_join(
    x = dbspec.adult %>% filter(forest_type=="evergreen") %>%
      group_by(treat,plot_ID) %>% summarise(nstem = n()),
    y = dbspec.adult %>% filter(forest_type=="evergreen", phenologyFIN=="D") %>%
      group_by(treat,plot_ID) %>% summarise(decistem = n())) %>% 
  mutate(deciprop = decistem/nstem) %>% 
  replace_na(list(deciprop = 0))

dbpsumm$decifrac_adult_evergreen <- 
  dbp$decifrac_adult_evergreen %>% group_by(treat) %>%
  summarise(y = mean(deciprop), ns = n(), ysd = sd(deciprop)) %>% 
  mutate(ymin = y - 1.96*(ysd/sqrt(ns)), ymax = y + 1.96*(ysd/sqrt(ns)),
         forest_type = "evergreen")


#### AGC = 0.5*AGB ~ logging frequency ####

dbp$AGB_deciduous <-
  dbagc.adult %>% 
   filter(forest_type=="deciduous") %>% group_by(treat,plot_ID) %>%
   summarise(agc = sum(agc_Mg)) %>% mutate(tagc = agc/(0.05))

dbpsumm$AGB_deciduous <-
  dbp$AGB_deciduous %>% mutate(tagc = log10(tagc)) %>% group_by(treat) %>% 
  summarise(y = mean(tagc), ns = n(), ysd = sd(tagc)) %>%
  mutate(ymin = y - 1.96*(ysd/sqrt(ns)), ymax = y + 1.96*(ysd/sqrt(ns)),
         forest_type = "deciduous")

dbp$AGB_evergreen <- 
  dbagc.adult %>%
  filter(forest_type=="evergreen") %>% group_by(treat,plot_ID) %>%
  summarise(agc = sum(agc_Mg)) %>% mutate(tagc = agc/(0.05))

dbpsumm$AGB_evergreen <-
  dbp$AGB_evergreen %>% mutate(tagc = log10(tagc)) %>% group_by(treat) %>% 
  summarise(y = mean(tagc), ns = n(), ysd = sd(tagc)) %>%
  mutate(ymin = y - 1.96*(ysd/sqrt(ns)), ymax = y + 1.96*(ysd/sqrt(ns)),
         forest_type = "evergreen")

#### composition difference - between forest type and logging treatment ####

bcdist <- 
  dcast(data = dbspec.adult %>% select(forest_type,plot_ID,code),
        formula = plot_ID + forest_type ~ code,
        fun.aggregate = length, value.var = "code") %>%
  column_to_rownames("plot_ID") %>% 
  select(-forest_type) %>% 
  decostand(x = ., method = "hellinger")

nmdsMod <- 
  metaMDS(comm = bcdist,
          distance = "bray",
          k = 2,
          try = 100, trymax = 300)

points.nmds <- 
  nmdsMod$points[,c(1:2)] %>% cbind() %>% as.data.frame()

compo.nmds <- 
  full_join(x = points.nmds %>% rownames_to_column("plot_ID"),
            y = db %>% select(plot_ID,treat,forest_type) %>% distinct(),
            by = "plot_ID")


#### running models & extracting coeffs (Table 1), p-value, F-statistics & df ####

#stem
modop[1:4] <- lapply(X = dbp[1:4],
                     FUN = function(x) 
                       glm(formula = nstem ~ treat,data = x,family = "poisson"))
modeldf[1:4] <- lapply(X = dbp[1:4],
                       FUN = function(x) 
                         glm(formula = nstem ~ treat,data = x,family = "poisson") %>% 
                         anova() %>% .$`Resid. Df` %>% .[2])
modelF[1:4] <- lapply(X = dbp[1:4],
                      FUN = function(x) 
                        glm(formula = nstem ~ treat,data = x,family = "poisson") %>% 
                        anova() %>% .$Deviance %>% .[2])

#canopy
modop[5:8] <- lapply(X = dbp[5:8],
                       FUN = function(x)
                         betareg(formula = ccover ~ treat,
                                 data = x %>% 
                                   mutate(ccover = ifelse(ccover==1,0.999,
                                                          ifelse(ccover==0,
                                                                 0.001,ccover)))))

#species richness
modop[9:12] <- lapply(X = dbp[9:12],
                     FUN = function(x)
                       glm(formula = round(expH) ~ treat,data = x,family = "poisson"))

modeldf[9:12] <- lapply(X = dbp[9:12],
                       FUN = function(x)
                         glm(formula = round(expH) ~ treat,data = x,family = "poisson") %>% 
                         anova() %>% .$`Resid. Df` %>% .[2])
modelF[9:12] <- lapply(X = dbp[9:12],
                      FUN = function(x) 
                        glm(formula = round(expH) ~ treat,data = x,family = "poisson") %>% 
                        anova() %>% .$Deviance %>% .[2])

#deciprop
modop[13:16] <- lapply(X = dbp[13:16],
                      FUN = function(x)
                        betareg(formula = deciprop ~ treat,data = x %>% 
                                           mutate(deciprop = ifelse(deciprop==1,0.999,
                                                                    ifelse(deciprop==0,0.001,deciprop)))))
                       
#agc
  modop[17:18] <- lapply(X = dbp[17:18],
                      FUN = function(x) glm(formula = log10(tagc) ~ treat,
                                                          data = x,family = "gaussian"))
modeldf[17:18] <- lapply(X = dbp[17:18],
                         FUN = function(x) anova(glm(formula = log10(tagc) ~ treat,
                                                     data = x,family = "gaussian"))$`Resid. Df`[2])
 modelF[17:18] <- lapply(X = dbp[17:18],
                        FUN = function(x) anova(glm(formula = log10(tagc) ~ treat,
                                                    data =  x,family = "gaussian"))$Deviance[2])

#compo
permadb <- full_join(y = bcdist %>% rownames_to_column("plot_ID"),
                     x = db %>% select(plot_ID,treat,forest_type) %>% distinct(),
                     by = "plot_ID")
rownames(permadb)<-permadb$plot_ID

permadb.contrasts <- permadb %>% select(plot_ID,treat,forest_type)

perma1 <- adonis(formula = permadb %>% select(-plot_ID,-treat,-forest_type) ~ treat,
                 data = permadb.contrasts,
                 method = "bray", permutations = 999)$aov.tab

perma2 <- adonis(formula = permadb %>% select(-plot_ID,-treat,-forest_type) ~ forest_type,
                 data = permadb.contrasts,
                 method = "bray", permutations = 999)$aov.tab

perma3 <- adonis(formula = permadb %>% select(-plot_ID,-treat,-forest_type) ~ forest_type + treat,
                 data = permadb.contrasts,
                 method = "bray", permutations = 999)$aov.tab

perma3 <- adonis(formula = permadb %>% select(-plot_ID,-treat,-forest_type) ~ forest_type*treat,
                 data = permadb.contrasts,
                 method = "bray", permutations = 999)$aov.tab

mop1<- jtools::export_summs(modop,digits.p = 2,
                            model.names = names(modop),
                            to.file = "docx",file.name = "October 2020 - Major revision/model_output.docx")

mop3<-data.frame(Fstat = round(do.call("rbind",modelF),4)) %>% rownames_to_column("modelname")
mop4<-data.frame(df = round(do.call("rbind",modeldf),4)) %>% rownames_to_column("modelname")


# calculating predicted % (deci-frac and canopy)
betareg_predicted <- 
modop[c(5:8,13:16)] %>% 
  map(.f = function(x) c(plogis(x$coefficients$mean[1]),
                         plogis(x$coefficients$mean[1] + x$coefficients$mean[2]),
                         plogis(x$coefficients$mean[1] + x$coefficients$mean[3])) %>%
        unname()) %>% 
  bind_rows()

betareg_perchange <- 
betareg_predicted %>% 
  add_row(100*(betareg_predicted[2,] - betareg_predicted[1,])) %>%  # L1 v/s B % change
  add_row(100*(betareg_predicted[3,] - betareg_predicted[1,]))      # L2 v/s B % change

# write.csv(x = betareg_perchange[-c(1:3),],
#           file = "October 2020 - Major revision/betareg_percentChange.csv",
#           row.names = F)   
# 
# write.csv(x = full_join(mop3,mop4),
#           file = "October 2020 - Major revision/model_output_Fdf.csv",
#           row.names = F) 
## table 1A - model outputs ##

# sink("October 2020 - Major revision/table1B_modelOP_permanova.csv")
# perma1 %>% round(.,3) %>%  as.data.frame()
# perma2 %>% round(.,3) %>%  as.data.frame()
# perma3 %>% round(.,3) %>%  as.data.frame()
# perma4 %>% round(.,3) %>%  as.data.frame()
# sink()


##### plotting figures 2, 3 and 4 ####

fig <- modop
labelsu <- c("B","L1","L2")

# dbpsumm <- 
# dbpsumm %>% 
#   map(.f = function(x) x %>% mutate(treat = case_when(treat=="B" & forest_type=="deciduous" ~ "B"#\n(12)",
#                                                       treat=="B" & forest_type=="evergreen" ~ "B",#\n(12)",
#                                                       treat=="L1" & forest_type=="deciduous" ~ "L1",#\n(18)",
#                                                       treat=="L1" & forest_type=="evergreen" ~ "L1",#\n(11)",
#                                                       treat=="L2" & forest_type=="deciduous" ~ "L2",#\n(12)",
#                                                       treat=="L2" & forest_type=="evergreen" ~ "L2"#\n(7)"
#                                                       ))
#       )
# 

fig[1:4] <-
  lapply(X = dbpsumm[1:4],
         FUN = 
           function(x)
           ggplot(data = x,
                  aes(y = y, x = treat)) +
           geom_pointrange(aes(ymin = ymin, ymax = ymax),
                           colour = "gray10", linetype = 2, shape = ifelse(x$forest_type=="deciduous",1,16)) +
           theme_light() + theme(plot.title = element_text(hjust = 0.5),axis.title = element_blank()) + 
           ggtitle("") + ylim(0,36))
           

tmpdat <- data.frame(x = c(1.875, 2.875), xend = c(2.125, 3.125), y = c(35,35),ann = c(".","."))
fig[[1]] <- fig[[1]] + geom_signif(data = tmpdat %>% mutate(ann = c("","*")), stat="identity", aes(x=x, xend=xend, y=y, yend=y, annotation=ann))
fig[[2]] <- fig[[2]] + geom_signif(data = tmpdat %>% mutate(ann = c("","")), stat="identity", aes(x=x, xend=xend, y=y, yend=y, annotation=ann))
fig[[3]] <- fig[[3]] + geom_signif(data = tmpdat %>% mutate(ann = c("***","")), stat="identity", aes(x=x, xend=xend, y=y, yend=y, annotation=ann))
fig[[4]] <- fig[[4]] + geom_signif(data = tmpdat %>% mutate(ann = c("","**")), stat="identity", aes(x=x, xend=xend, y=y, yend=y, annotation=ann))

fig[5:8] <-
  lapply(X = dbpsumm[5:8],
         FUN = 
           function(x)
             ggplot(data = x, aes(y = y, x = treat)) +
           geom_pointrange(aes(ymin = ymin, ymax = ymax),
                           colour = "gray10", linetype = 2, shape = ifelse(x$forest_type=="deciduous",1,16)) +
           theme_light() + theme(plot.title = element_text(hjust = 0.5),axis.title = element_blank()) + 
           ggtitle("") + ylim(0,1))

tmpdat <- data.frame(x = c(1.875, 2.875), xend = c(2.125, 3.125), y = c(0.9,0.9),ann = c(".","."))
fig[[5]] <- fig[[5]] + geom_signif(data = tmpdat %>% mutate(ann = c("***","***")),stat="identity", aes(x=x, xend=xend, y=y, yend=y, annotation=ann))
fig[[6]] <- fig[[6]] + geom_signif(data = tmpdat %>% mutate(ann = c("***","***")),stat="identity", aes(x=x, xend=xend, y=y, yend=y, annotation=ann))
fig[[7]] <- fig[[7]] + geom_signif(data = tmpdat %>% mutate(ann = c("","**")), stat="identity", aes(x=x, xend=xend, y=y, yend=y, annotation=ann))
fig[[8]] <- fig[[8]] + geom_signif(data = tmpdat %>% mutate(ann = c("",".")), stat="identity", aes(x=x, xend=xend, y=y, yend=y, annotation=ann))

fig[9:12] <-
  lapply(X = dbpsumm[9:12],
         FUN = 
           function(x)
             ggplot(data = x, aes(y = y, x = treat)) +
           geom_pointrange(aes(ymin = ymin, ymax = ymax),
                           colour = "gray10", linetype = 2, shape = ifelse(x$forest_type=="deciduous",1,16)) +
           theme_light() + theme(plot.title = element_text(hjust = 0.5),axis.title = element_blank()) + 
           ggtitle("") + ylim(0,20))

tmpdat <- data.frame(x = c(1.875, 2.875), xend = c(2.125, 3.125), y = c(18,18),ann = c(".","."))
fig[[9]] <- fig[[9]] + geom_signif(data = tmpdat %>% mutate(ann = c("***",".")),stat="identity", aes(x=x, xend=xend, y=y, yend=y, annotation=ann))
fig[[10]] <- fig[[10]] + geom_signif(data = tmpdat %>% mutate(ann = c("","")),stat="identity", aes(x=x, xend=xend, y=y, yend=y, annotation=ann))
fig[[11]] <- fig[[11]] + geom_signif(data = tmpdat %>% mutate(ann = c("","**")), stat="identity", aes(x=x, xend=xend, y=y, yend=y, annotation=ann))
fig[[12]] <- fig[[12]] + geom_signif(data = tmpdat %>% mutate(ann = c("",".")), stat="identity", aes(x=x, xend=xend, y=y, yend=y, annotation=ann))

fig[13:16] <-
  lapply(X = dbpsumm[13:16],
         FUN = 
           function(x)
             ggplot(data = x, aes(y = y, x = treat)) +
           geom_pointrange(aes(ymin = ymin, ymax = ymax),
                           colour = "gray10", linetype = 2, shape = ifelse(x$forest_type=="deciduous",1,16)) +
           theme_light() + theme(plot.title = element_text(hjust = 0.5),axis.title = element_blank()) + 
           ggtitle("") + ylim(-0.02,0.62))

tmpdat <- data.frame(x = c(1.875, 2.875), xend = c(2.125, 3.125), y = c(0.6,0.6),ann = c(".","."))
fig[[13]] <- fig[[13]] + geom_signif(data = tmpdat %>% mutate(ann = c("","")), stat="identity", aes(x=x, xend=xend, y=y, yend=y, annotation=ann))
fig[[14]] <- fig[[14]] + geom_signif(data = tmpdat %>% mutate(ann = c("","")), stat="identity", aes(x=x, xend=xend, y=y, yend=y, annotation=ann))
fig[[15]] <- fig[[15]] + geom_signif(data = tmpdat %>% mutate(ann = c("**",".")), stat="identity", aes(x=x, xend=xend, y=y, yend=y, annotation=ann))
fig[[16]] <- fig[[16]] + geom_signif(data = tmpdat %>% mutate(ann = c(".","")), stat="identity", aes(x=x, xend=xend, y=y, yend=y, annotation=ann))

fig[17:18] <-
  lapply(X = dbpsumm[17:18],
         FUN = 
           function(x)
             ggplot(data = x %>% 
                      mutate(y = 10^y, ysd = 10^ysd, ymin = 10^ymin, ymax = 10^ymax),
                    aes(y = y, x = treat)) +
           geom_pointrange(aes(ymin = ymin, ymax = ymax),
                           colour = "gray10", linetype = 2, shape = ifelse(x$forest_type=="deciduous",1,16)) +
           theme_light() + theme(plot.title = element_text(hjust = 0.5),axis.title = element_blank()) + 
           ggtitle("") + ylim(0,280))

tmpdat <- data.frame(x = c(1.875, 2.875), xend = c(2.125, 3.125), y = c(260,260),ann = c(".","."))
fig[[17]] <- fig[[17]] + geom_signif(data = tmpdat %>% mutate(ann = c("","***")), stat="identity", aes(x=x, xend=xend, y=y, yend=y, annotation=ann))
fig[[18]] <- fig[[18]] + geom_signif(data = tmpdat %>% mutate(ann = c("","")), stat="identity", aes(x=x, xend=xend, y=y, yend=y, annotation=ann))

fig[[19]] <-
  ggplot(data = compo.nmds,
         aes(y=MDS2,x=MDS1)) + 
  geom_point(aes(colour = forest_type,shape = treat),
             size = 3) + 
  coord_equal() +
  theme_light() + theme(legend.position = "top",
                        plot.title = element_text(hjust = 0.5)) +
  scale_colour_manual(name = "", values = c("dodgerblue2","forestgreen")) + 
  scale_shape_manual(name = "", values = c(15,0,7)) +
  annotate("text", x=-1.75, y=1.2, size=4, hjust=0,
           label = expression("composition ~ logging treatment** R"^{2}*"=0.068")) +
  annotate("text", x=-1.75, y=1.1, size=4, hjust=0,
           label = expression("composition ~ forest type*** R"^{2}*"=0.212"))

fig2_toprow<-
  cowplot::plot_grid(plotlist = fig[c("stem_pole_deciduous","stem_pole_evergreen",
                                      "stem_adult_deciduous","stem_adult_evergreen")],
                       nrow = 1,ncol = 4,
                       labels = c("pole-tree density","",
                                  "adult-tree density",""),label_fontface = "plain",label_size = 10)

fig2_botrow<-
  cowplot::plot_grid(plotlist = fig[c("understory_deciduous","understory_evergreen",
                                        "overstory_deciduous","overstory_evergreen",
                                        "AGB_deciduous","AGB_evergreen")],
                       nrow = 1,ncol = 6, 
                       labels = c("understory cover","","overstory cover","","above-ground carbon",""),label_fontface = "plain",label_size = 10)


fig3<-
  cowplot::plot_grid(plotlist = fig[c("richness_pole_deciduous","richness_pole_evergreen",
                                      "richness_adult_deciduous","richness_adult_evergreen",
                                      "decifrac_pole_deciduous","decifrac_pole_evergreen",
                                      "decifrac_adult_deciduous","decifrac_adult_evergreen")],
                     nrow = 2, ncol = 4,
                labels = c("pole-tree richness","","adult-tree richness","","pole-tree deci-fraction","","adult-tree deci-fraction"),label_fontface = "plain",label_size = 10)
fig4<-fig[[19]]


for(i in 1:18)
{
  ggsave(plot = fig[[i]],filename = paste0("figures/individual_fig",i,"_",names(fig)[i],".png"),
         dpi = 600, units = "in", height = 3, width = 1.2)
}

legendplot <- 
  get_legend(plot = dbpsumm[1:2] %>% bind_rows() %>% ggplot(data = .) + geom_point(aes(y=y,x=treat,shape=forest_type)) + scale_shape_manual(name = "", values = c(1,16)) + theme_light())

ggsave(plot = legendplot,
       filename = "figures/legend.png",
       device = "png", width = 4,height = 2.5,units = "in", dpi = 600)

 ggsave(plot = fig2_toprow,
        filename = "figures/fig2_toprow_stemden.png",
        device = "png", width = 4,height = 2.5,units = "in", dpi = 600)
 ggsave(plot = fig2_botrow,
        filename = "figures/fig2_bottomrow_canopy_agb.png",
        device = "png", width = 6,height = 2.5,units = "in", dpi = 600)
 ggsave(plot = fig3,
        filename = "figures/fig3_specrich_decifrac.png",
        device = "png", width = 4,height = 5,units = "in", dpi = 600)
 ggsave(plot = fig4,
        filename = "figures/fig4_compo.png",
        device = "png", width = 7, height = 7,units = "in", dpi = 600)

#### SUPPLEMENTARY FIGURES and tables ######
# see AndamansMScAS2020_supplementary.R

