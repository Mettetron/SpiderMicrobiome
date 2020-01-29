library(ggplot2)
library(dplyr)
setwd("~//")
# import data and trim white space off character columns
spider.data <- read.csv("spider_symbionts.csv", stringsAsFactors = FALSE)
truth <- sapply(spider.data, is.character)
spider.data <- data.frame(cbind(sapply(spider.data[, truth], trimws, which="both"), spider.data[,!truth]))

# sum data by spider species _ symbiont
spider.data$species_bug <- paste0(spider.data$species, "_", spider.data$Symbiont)
spider.data$incidence_count <- round(as.numeric(as.character(spider.data$incidence))*as.numeric(as.character(spider.data$n)), 0)

spider.dat.species.summed <- spider.data %>% 
  group_by(species_bug) %>% 
  summarise(incidence_count=sum(incidence_count), n=sum(n))

s.dat.sum <- as.data.frame(spider.dat.species.summed)
s.dat.sum$studies <- table(spider.data$species_bug)
table(spider.data$species_bug)
table(spider.dat.species.summed$species_bug)

s.dat.sum$incidence_frac <- round(s.dat.sum$incidence_count/s.dat.sum$n, 2)
hej <- strsplit(s.dat.sum$species_bug, "_")
s.dat.sum$species <- sapply(strsplit(s.dat.sum$species_bug, "_"), `[`, 1)
s.dat.sum$symbiont <- sapply(strsplit(s.dat.sum$species_bug, "_"), `[`, 2)
s.dat.sum$genus <- sapply(strsplit(s.dat.sum$species, " "), `[`, 1)
s.dat.sum$family <- spider.data$family[match(s.dat.sum$species, spider.data$species)]

# get total number of individuals sampled for each study
ind_spec <- spider.data[, c("ref", "species", "n")]
ind_spec <- ind_spec[!duplicated(ind_spec), ]
ind_spec_summed <- ind_spec %>%
  group_by(species) %>%
  summarise(n=sum(n))

s.dat.sum$ind_per_species <- ind_spec_summed$n[match(s.dat.sum$species, ind_spec_summed$species)]
s.dat.sum$species <- paste0(s.dat.sum$species, " (", s.dat.sum$ind_per_species, ")")

# put the families in rough phylogenetic order
family.order <- c("Eresidae", "Theridiidae", "Tetragnathidae", "Linyphiidae", 
                  "Nephilidae", "Araneidae", "Oecobiidae", "Uloboridae", 
                  "Titanoecidae", "Urocteidae", "Deinopidae", "Amaurobiidae", "Cybaeidae",  
                  "Agelenidae","Desidae", "Philodromidae","Miturgidae", 
                  "Salticidae", "Thomisidae", "Clubionidae", "Gnaphosidae", 
                  "Oxyopidae", "Lycosidae", "Pisauridae", "Dysderidae",
                  "Pholcidae", "Telemidae", "Scytodidae")

#unique(s.dat.sum$family)[!unique(s.dat.sum$family) %in% family.order]
s.dat.sum$family <- factor(s.dat.sum$family, levels=family.order)
s.dat.sum$fam.num <- as.numeric(s.dat.sum$family)

s.dat.sum <- s.dat.sum[order(s.dat.sum$fam.num), ]

s.dat.sum$species <- factor(s.dat.sum$species, levels = unique(s.dat.sum$species))

s.dat.sum$symbiont <- factor(s.dat.sum$symbiont, levels=c("Wolbachia", "Rickettsia", "Cardinium", "Spiroplasma"))


p <- ggplot(s.dat.sum, aes(as.factor(symbiont), species, group=species)) +
  theme_bw() + 
  geom_tile(aes(fill = incidence_frac)) + 
  geom_text(aes(label = paste(incidence_count, "/", n)), size=2) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1.1, vjust=1.1, size = 10, color = 'black')) +
  scale_fill_gradient(low = "blue", high = "red", na.value = "black", name="Prevalence\nin species") +
  facet_grid(family ~ ., scales="free", space="free") +
  theme(strip.text.y = element_text(size = 10, angle = 0)) +
  theme(panel.spacing=unit(0.1, "lines")) + # change spacing between facets
  xlab("Symbiont") +
  theme(legend.position="bottom")
  
ggsave("endosymbHeatmap.png", p,
       scale = 1, width = 12, height = 100, units = "cm",
       dpi = 300, limitsize = TRUE)

# calculate amount of species infected by each symbionts
symb.vect <- c()
species.count <- c()
for (symb in unique(s.dat.sum$symbiont)) {
  #symb <- "Cardinium"
  df.symb <- s.dat.sum[s.dat.sum$symbiont == symb, ]
  symb.vect <- c(symb.vect, sum(df.symb$incidence_count >= 1))
  species.count <- c(species.count, nrow(df.symb))
}

results.df <- data.frame(symbiont = unique(s.dat.sum$symbiont), 
                         species.with.symb = symb.vect, 
                         species.total = species.count)


# Genus summed heatmap ----------------------------------------------------

spider.data$genus <- sapply(strsplit(as.character(spider.data$species), " "), `[`, 1)
spider.data$genus_bug <- paste0(spider.data$family,"_", spider.data$genus, "_", spider.data$Symbiont)

s.dat.sum.genus <- spider.data %>% 
  group_by(genus_bug) %>% 
  summarise(incidence_count=sum(incidence_count), n=sum(n))

s.dat.sum.genus <- as.data.frame(s.dat.sum.genus)

s.dat.sum.genus$incidence_frac <- round(s.dat.sum.genus$incidence_count/s.dat.sum.genus$n, 2)
s.dat.sum.genus$family <- sapply(strsplit(s.dat.sum.genus$genus_bug, "_"), `[`, 1)
s.dat.sum.genus$genus <- sapply(strsplit(s.dat.sum.genus$genus_bug, "_"), `[`, 2)
s.dat.sum.genus$symbiont <- sapply(strsplit(s.dat.sum.genus$genus_bug, "_"), `[`, 3)

# find number of species in each genus
spec.gen <- spider.data[, c("genus", "species")]
spec.gen <- spec.gen[!duplicated(spec.gen), ]
#table(spec.gen$genus)
s.dat.sum.genus$num_spec <- table(spec.gen$genus)[s.dat.sum.genus$genus]


s.dat.sum.genus$family <- factor(s.dat.sum.genus$family, levels=family.order)
s.dat.sum.genus$symbiont <- factor(s.dat.sum.genus$symbiont, levels=c("Wolbachia", "Rickettsia", "Cardinium", "Spiroplasma"))

s.dat.sum.genus$genus <- paste0(s.dat.sum.genus$genus, " (", s.dat.sum.genus$num_spec, ")")


q <- ggplot(s.dat.sum.genus, aes(as.factor(symbiont), genus, group=genus)) +
  theme_bw() + 
  geom_tile(aes(fill = incidence_frac)) + 
  geom_text(aes(label = paste(incidence_count, "/", n)), size=2) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1.1, vjust=1.1, size = 15, color = 'black')) +
  scale_fill_gradient(low = "blue", high = "red", na.value = "black", name="Prevalence\nin genus") +
  facet_grid(family ~ ., scales="free", space="free") +
  theme(strip.text.y = element_text(size = 10, angle = 0)) +
  theme(panel.spacing=unit(0.1, "lines")) + # change spacing between facets
  theme(legend.position="bottom") +
  ylab("Spider genus") +
  theme(axis.title.x = element_blank()) 

ggsave("endosymbHeatmap_GENUS.pdf", q,
       scale = 1, width = 12, height = 50, units = "cm",
       dpi = 300, limitsize = TRUE)


##### split plot

plot1.fam <- c("Eresidae", "Theridiidae", "Tetragnathidae", "Linyphiidae", 
               "Nephilidae", "Oecobiidae", "Uloboridae","Titanoecidae")

plot2.fam <- c( "Araneidae",
                "Urocteidae", "Deinopidae", "Amaurobiidae", "Cybaeidae",  
                "Agelenidae","Desidae", "Philodromidae","Miturgidae", 
                "Salticidae", "Thomisidae", "Clubionidae", "Gnaphosidae", 
                "Oxyopidae", "Lycosidae", "Pisauridae", "Dysderidae",
                "Pholcidae", "Telemidae", "Scytodidae")

library(scales)
library(ggplot2)
dat1 <- s.dat.sum.genus[s.dat.sum.genus$family %in% plot1.fam, ]
plot1 <- ggplot(dat1, aes(as.factor(symbiont), genus, group=genus)) +
  theme_bw() + 
  geom_tile(aes(fill = incidence_frac)) + 
  geom_text(aes(label = paste(incidence_count, "/", n)), size=2) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1.1, vjust=1.1, size = 10, color = 'black')) +
  #scale_fill_gradient(low = "blue", high = "red", labels = percent, na.value = "black", name="Prevalence\nin genus") +
  scale_fill_gradient(low = "blue", high = "red", na.value = "black", name="Prevalence\nin genus") +
  facet_grid(family ~ ., scales="free", space="free") +
  theme(strip.text.y = element_text(size = 10, angle = 0)) +
  theme(panel.spacing=unit(0.1, "lines")) + # change spacing between facets
  # theme(text=element_text(family="Helvetica"), axis.title.y = element_blank(), axis.title.x=element_blank()) +
  theme(legend.position="bottom") +
  #ylab("Spider genus") +
  theme(axis.title = element_blank()) 

# ggsave("plotforlegend.pdf", plot1,
#        scale = 1, width = 10, height = 25, units = "cm",
#        dpi = 300, limitsize = TRUE)

dat2 <- s.dat.sum.genus[s.dat.sum.genus$family %in% plot2.fam, ]
plot2 <- ggplot(dat2, aes(as.factor(symbiont), genus, group=genus)) +
  theme_bw() + 
  geom_tile(aes(fill = incidence_frac)) + 
  geom_text(aes(label = paste(incidence_count, "/", n)), size=2) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1.1, vjust=1.1, size = 10, color = 'black')) +
  scale_fill_gradient(low = "blue", high = "red", na.value = "black", name="Prevalence\nin genus") +
  facet_grid(family ~ ., scales="free", space="free") +
  theme(strip.text.y = element_text(size = 10, angle = 0)) +
  theme(panel.spacing=unit(0.1, "lines")) + # change spacing between facets
  # theme(text=element_text(family="Helvetica"), axis.title.y = element_blank(), axis.title.x=element_blank()) +
  theme(legend.position="none") +
  #ylab("Spider genus") +
  theme(axis.title = element_blank()) 


library(gridExtra)

p3 <- grid.arrange(rbind(ggplotGrob(plot1), ggplotGrob(plot2)))
                   
ggsave("endosymbHeatmap_GENUS_2.pdf", p3,
       scale = 1, width = 10, height = 50, units = "cm",
       dpi = 300, limitsize = TRUE)




# make table for thesis ---------------------------------------------------


substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

spider.data$study.year <- substrRight(as.character(spider.data$ref), 5)
spider.data$study.year <- substr(spider.data$study.year, 1, 4)

species.vector <- c()
species.num <- c()
symbiont.vector <- c()
for (study in unique(spider.data$ref)){
  species <- sort(unique(spider.data[spider.data$ref == study, "species"]))
  species.num <- c(species.num, length(species))
  species.vector <- c(species.vector, paste(species, collapse=". "))
  symb <- unique(spider.data[spider.data$ref == study, "Symbiont"])
  symb <- unique(spider.data[spider.data$ref == study, "Symbiont"])[order(match(unique(spider.data[spider.data$ref == study, "Symbiont"]), c("Wolbachia", "Rickettsia", "Cardinium", "Spiroplasma")))]
  symbiont.vector <- c(symbiont.vector, paste(symb, collapse=". "))
}

study.sum.df <- spider.data[, c("ref", "method", "study.year")]
study.sum.df <- study.sum.df[!duplicated(study.sum.df), ]



# get n
spider.data$study_species_n <- paste0(spider.data$ref, "_", spider.data$species, "_", spider.data$n)

n.df <- spider.data %>%
  group_by(study_species_n) %>%
  summarise(n.cheat=sum(n))


n.df <- as.data.frame(n.df)
n.df$study <- sapply(strsplit(n.df$study_species_n, "_"), `[`, 1)
n.df$species <- sapply(strsplit(n.df$study_species_n, "_"), `[`, 2)
n.df$n <- sapply(strsplit(n.df$study_species_n, "_"), `[`, 3)



n.df.final <- n.df %>%
  group_by(study) %>%
  summarise(n.study = sum(as.numeric(n)))

# put n in
study.sum.df$n <- n.df.final$n.study[match(as.character(study.sum.df$ref), n.df.final$study)]

if (identical(study.sum.df$ref, unique(spider.data$ref))) {
  study.sum.df$species <- species.vector
  study.sum.df$numberOfSpecies <- species.num
  study.sum.df$symbiont <- symbiont.vector
  
}

study.sum.df <- study.sum.df[order(study.sum.df$study.year), ]


write.csv(study.sum.df, "SpiderSymbMeta_studies.csv")
# p3 <- grid.arrange(cbind(ggplotGrob(dum.p + theme(plot.margin = unit(c(5.1, 2, 4.1, 5.1), "pt")) +
#                                       scale_x_continuous(breaks=seq(0, 30, 5))), 
#                          ggplotGrob(mim.p + theme(plot.margin = unit(c(5.1, 2, 4.1, 1), "pt"), 
#                                                   axis.text.y = element_blank(), 
#                                                   axis.ticks.y = element_blank())), 
#                          ggplotGrob(sar.p + theme(plot.margin = unit(c(5.1, 4.1, 4.1, 1), "pt"), 
#                                                   axis.text.y = element_blank(), 
#                                                   axis.ticks.y = element_blank()))))
#   