# Impact of repeat logging on forest recovery
# Akshay Surendra, Anand M Osuri, Jayashree Ratnam
# August 2020

#### Packages ######

# rm(list = ls())
library(vegan)
library(tidyverse)
library(gridExtra)
library(reshape2)
library(grid)
library(ggsignif)
library(cowplot)
library(betareg)
library(extrafont)
# font_import()
# loadfonts(device="win")
fonts()

### Table 1 #####

db_fgeo <- read_csv("plot_characteristics.csv")

tab1 <-
  full_join(x = db_fgeo %>% group_by(LoggingTreatment,forest_type) %>% summarise(Nplots = n_distinct(plotID)),
            y = db_fgeo %>%
              group_by(LoggingTreatment,forest_type) %>%
              summarise_at(.vars = c("temperatureBIOCLIM_degC","rainfallBIOCLIM_mmPerYr","elevSRTM30_m","slopeSRTM30_degree","distSettlement_m"),
                           .funs = function(x) paste0(round(median(x,na.rm=T),1),
                                                      " (",round(min(x,na.rm=T),1),
                                                      "-",round(max(x,na.rm=T),1),")")),
            by = c("forest_type","LoggingTreatment"))

write_csv(x = tab1,path = "output/Table1_plotchar.csv")


#### Data input and clean-up ####

db0 <- read_csv("plotdata.csv")
phenology <- read_csv("treecodes.csv")
canopy_raw <- read_csv("canopy.csv")

removeplots <- c("BANL0303","BAT1P10","MAT1C6P09","MAT2C2P09")
### Identifying and removing plots that are erroneously classified, (see supplementary.R script)

dbtmp1 <-
  db0 %>%
  filter(!plot_ID %in% removeplots) %>%
  mutate(treat = case_when(treatment ==      "TL" ~ "L2",
                           treatment == "OL_near" ~ "L1",
                           treatment ==  "OL_far" ~ "B",
                           treatment =="baseline" ~ "B"),
         treat = factor(treat, levels = c("B","L1","L2")),
         forest_type = factor(forest_type, levels = c("deciduous","evergreen")),
         plot_ID = factor(plot_ID)) %>%
  select(-treatment) %>%
  mutate(agb = data.frame(agbstem1 = ( (0.0673) * ((dbtmp$wd_m) * (dbtmp$dia1^2) * (dbtmp$height_m))^(0.976) ),
                          agbstem2 = ( (0.0673) * ((dbtmp$wd_m) * (dbtmp$dia2^2) * (dbtmp$height_m))^(0.976) ),
                          agbstem3 = ( (0.0673) * ((dbtmp$wd_m) * (dbtmp$dia3^2) * (dbtmp$height_m))^(0.976) )) %>%
           apply(X = ., MARGIN = 1, FUN = sum, na.rm = T),
         agc_Mg = 0.5*agb/1000) %>% #Chave 2014 pan-tropical formula used; Carbon fraction ~50% of AGB (Chave 2005)
  left_join(x = .,
            y = phenology %>% select(species_ID,phenologyFIN))

setdiff(x = unique(db$species_ID),
        y = unique(phenology$species_ID)) # only stump has no phenology

cc <-
  canopy_raw %>%
  mutate(lowerstorey = apply(X = cc[,c(2:10)],
                             MARGIN = 1,FUN = sum))
cc$lowerstorey <- apply(X = cc[,c(2:10)],
                         MARGIN = 1,FUN = sum)
cc$upperstorey <- apply(X = cc[,c(11:19)],
                         MARGIN = 1,FUN = sum) ## total canopy over 45 metres summed up by adding counts from 9 5-metre sections

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
  mutate(ccover = ccover/100) # percentage to proportion

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

#### canopy cover ~ logging frequency (understory and overstory) ####

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
                            to.file = "docx",file.name = "output/Table_2_temp.docx")

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

 write.csv(x = betareg_perchange[-c(1:3),],
           file = "output/Table2_Betareg_percentChange.csv",
           row.names = F)

 write.csv(x = full_join(mop3,mop4),
           file = "output/Table2_Fdf.csv",
           row.names = F)

 sink("output/tableS2_modelOP_permanova.csv")
 perma1 %>% round(.,3) %>%  as.data.frame()
 perma2 %>% round(.,3) %>%  as.data.frame()
 perma3 %>% round(.,3) %>%  as.data.frame()
 perma4 %>% round(.,3) %>%  as.data.frame()
 sink()

##### plotting figures 2, 3 and 4 ####

fig <- modop
labelsu <- c("B","L1","L2")

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

legendplot <- get_legend(plot = dbpsumm[1:2] %>% bind_rows() %>% ggplot(data = .) + geom_point(aes(y=y,x=treat,shape=forest_type)) + scale_shape_manual(name = "", values = c(1,16)) + theme_light())

###

title1 <- ggdraw() + draw_label("a) Pole-tree density \n (0.01ha)",fontface = "plain",fontfamily = "Times New Roman",size = 10)
title2 <- ggdraw() + draw_label("b) Adult-tree density \n (0.05ha)",fontface = "plain",fontfamily = "Times New Roman",size = 10)
title3 <- ggdraw() + draw_label("c) Understory canopy \n (0: no cover - 1: full cover)",fontface = "plain",fontfamily = "Times New Roman",size = 10)
title4 <- ggdraw() + draw_label("d) Overstory canopy \n (0: no cover - 1: full cover)",fontface = "plain",fontfamily = "Times New Roman",size = 10)
title5 <- ggdraw() + draw_label("e) Above-ground carbon \n (MgC per ha)",fontface = "plain",fontfamily = "Times New Roman",size = 10)

fig2 <-
plot_grid(
  plot_grid( # row 1
    ggplot() + theme_minimal(),
    plot_grid(title1,plot_grid(fig$stem_pole_deciduous,fig$stem_pole_evergreen), ncol = 1, rel_heights = c(0.15,0.85)),
    plot_grid(title2,plot_grid(fig$stem_adult_deciduous,fig$stem_adult_evergreen),ncol = 1, rel_heights = c(0.15,0.85)),
    ggplot() + theme_minimal(),
    ncol = 4,nrow=1, rel_widths = c(0.167,0.333,0.333,0.167)),
  ggdraw(legendplot), # intermediate row for the legend
  plot_grid( # row 2
    plot_grid(title3,plot_grid(fig$understory_deciduous,fig$understory_evergreen), ncol = 1, rel_heights = c(0.15,0.85)),
    plot_grid(title4,plot_grid(fig$overstory_deciduous,fig$overstory_evergreen), ncol = 1, rel_heights = c(0.15,0.85)),
    plot_grid(title5,plot_grid(fig$AGB_deciduous,fig$AGB_evergreen), ncol = 1, rel_heights = c(0.15,0.85)),
    ncol = 3,nrow=1,rel_widths = c(0.333,0.333,0.333)),
  nrow = 3, ncol = 1,rel_heights = c(0.42,0.16,0.42))

ggsave(filename = "output/figure2_andaman-logging.png",
       plot = fig2,device = "png",width = 8,height = 6,units = "in",dpi = 600)

#

title6 <- ggdraw() + draw_label("a) Pole-tree richness \n (0.01ha)",fontface = "plain",fontfamily = "Times New Roman",size = 10)
title7 <- ggdraw() + draw_label("b) Adult-tree richness \n (0.05ha)",fontface = "plain",fontfamily = "Times New Roman",size = 10)
title8 <- ggdraw() + draw_label("c) Pole-tree deciduous fraction \n (0.01ha)",fontface = "plain",fontfamily = "Times New Roman",size = 10)
title9 <- ggdraw() + draw_label("d) Adult-tree deciduous fraction \n (0.05ha)",fontface = "plain",fontfamily = "Times New Roman",size = 10)

fig3 <-
  plot_grid(
    plot_grid( # row 1
      plot_grid(title6,plot_grid(fig$richness_pole_deciduous,fig$richness_pole_evergreen), ncol = 1, rel_heights = c(0.15,0.85)),
      plot_grid(title7,plot_grid(fig$richness_adult_deciduous,fig$richness_adult_evergreen),ncol = 1, rel_heights = c(0.15,0.85)),
      ncol = 2,nrow=1, rel_widths = c(0.333,0.333)),
    ggdraw(plot = legendplot,xlim = c(0,1),ylim = c(0,1)), # intermediate row for the legend
    plot_grid( # row 2
      plot_grid(title8,plot_grid(fig$decifrac_pole_deciduous,fig$decifrac_pole_evergreen), ncol = 1, rel_heights = c(0.15,0.85)),
      plot_grid(title9,plot_grid(fig$decifrac_adult_deciduous,fig$decifrac_adult_evergreen), ncol = 1, rel_heights = c(0.15,0.85)),
      ncol = 2,nrow=1,rel_widths = c(0.333,0.333)),
    nrow = 3, ncol = 1,rel_heights = c(0.42,0.16,0.42))

ggsave(filename = "output/figure3_andaman-logging.png",
       plot = fig3,device = "png",width = 8,height = 6,units = "in",dpi = 600)

fig4<-fig[[19]]

ggsave(filename = "output/figure4_andaman-logging.png",
       plot = fig4,device = "png",width = 5,height = 5,units = "in",dpi = 600)


#### SUPPLEMENTARY FIGURES and tables - see AndamansMScAS2020_supplementary.R





