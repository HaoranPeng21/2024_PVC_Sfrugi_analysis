
##### 16S pipeline ####
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
devtools::install_github("jbisanz/qiime2R")
devtools::install_github("vmikk/metagMisc")
BiocManager::install("MicrobiotaProcess")
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(qiime2R)
library(metagMisc)
library(ggplotify)
library(cowplot)
library(vegan)
library(MicrobiotaProcess)
setwd("T.molitor/1.1--Qiime_4")

THMAC.ps <-qza_to_phyloseq(
  features="table-with-phyla-no-metochondria-no-chloroplast.qza",
  tree="rooted-tree.qza",
  taxonomy="taxonomy.qza",
  metadata = "metadata.txt")

THMAC.ps
THMAC.ps_filter <- subset_samples(THMAC.ps, Sites != "Cornfield soil")
THMAC.ps.domint <- filter_taxa(THMAC.ps_filter, function(x) sum(x>2) > 1,TRUE)
THMAC.ps.domint
THMAC.ps.domint <- subset_samples(THMAC.ps.domint, Sites == "Cornfield soil" )
THMAC.ps.domint
THMAC.ps.domint <- subset_samples(THMAC.ps.domint, Sites == "Intestinal faces" )
THMAC.ps.domint
THMAC.ps.domint <- subset_samples(THMAC.ps.domint, Sites == "Excreted faeces")

THMAC.ps.domint <- filter_taxa(THMAC.ps, function(x) sum(x>2) > 1,TRUE)
THMAC.ps.domint <- subset_samples(THMAC.ps.domint, Sites == "Intestinal mucosa")
THMAC.ps.domint

ps_sites <- THMAC.ps.domint %>% as.MPSE() 

ps_sites %<>%
  mp_cal_abundance( # for each samples
    .abundance = RareAbundance
  ) %>%
  mp_cal_abundance( # for each groups 
    .abundance=RareAbundance,
    .group=Places
  )
ps_sites
phylum1_colors <- c("p__Actinobacteriota" = "#9FBA95", "p__Bacteroidota" = "#E6CECF", "p__Firmicutes" = "#B696B6"
                    , "p__Proteobacteria" = "#80C1C4","p__Chloroflexi" = "#fae69e","p__Dadabacteria" = "#f2b56f",
                    "p__Fusobacteria" = "#BB9393","p__Patescibacteria" = "#808080","p__Spirochaetes" = "#7CAEF0",
                    "p__Tenericutes" = "#cee9dc","p__Verrucomicrobia" = "#ffadbb")
c("#9FBA95", "#E6CECF","#B696B6", "#80C1C4")
compo3.1 <- ps_sites %>%
  mp_plot_abundance(
    .abundance=RareAbundance,
    .group=Sites, 
    taxa.class = Genus, 
    topn = 10,
    relative = TRUE,
    plot.group = F
  ) +
  scale_fill_manual(
    values = c("#808080","lightgrey", "#ffadbb", "#cee9dc","#f2b56f", "#7CAEF0",
               "#fae69e", "#80C1C4", "#E6CECF", "#B696B6", "#9FBA95"),
    guide = guide_legend(keywidth = 0.6, keyheight = 0.6, label.theme = element_text(size = 10))
  ) +
  theme_classic() + 
  theme(panel.grid =element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 10,colour = "black"),
        legend.position = 'right')
compo3.1

sphingo_ef <- c(61.402411481,70.238715730,69.116369575,33.058827617)
mean(sphingo_ef)
sd(sphingo_ef)

ps_sites %<>% 
  mp_decostand(.abundance=Abundance)
ps_sites

ps_sites %<>% mp_cal_dist(.abundance=hellinger, distmethod="bray")
ps_sites

ps_sites %<>% 
  mp_cal_pcoa(.abundance=hellinger, distmethod="bray")
# The dimensions of ordination analysis will be added the colData slot (default).
ps_sites

ps_sites %<>%
  mp_adonis(.abundance=hellinger, .formula=~Places, distmethod="bray", permutations=9999, action="add")
ps_sites %>% mp_extract_internal_attr(name=adonis)

library(ggplot2)
### Alpha diversity
ps_sites %<>% 
  mp_cal_alpha(.abundance=RareAbundance)
ps_sites
alpha4 <- ps_sites %>% 
  mp_plot_alpha(
    .group=Places, 
    test = "wilcox.test",
    .alpha=c(Chao1, Shannon)
  ) +
  scale_fill_manual(values=c("#9FBA95", "#E6CECF","#B696B6", "#80C1C4"), guide="none") +
  scale_color_manual(values=c("#9FBA95", "#E6CECF","#B696B6", "#80C1C4"), guide="none")

#ps_sites[["Shannon"]]

alpha4


beta1 <- ps_sites %>% 
  mp_plot_ord(
    .ord = pcoa, 
    .group = Sites, 
     show.adonis = T,
    .color = Sites, 
    .size = Shannon, 
    .alpha = Shannon,
    ellipse = TRUE,
    show.legend = FALSE # don't display the legend of stat_ellipse 
  ) +
  scale_fill_manual(
    values = c("#9FBA95", "#E6CECF","#B696B6", "#80C1C4"), 
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
  ) +
  scale_color_manual(
    values=c("#9FBA95", "#E6CECF","#B696B6", "#80C1C4"),
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
  ) +
  scale_size_continuous(
    range=c(0.5, 3),
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
  ) #+
  #facet_wrap(~ Sites)
beta1

aplot::plot_list(gglist=list(compo1,compo2, compo3,compo4), tag_levels="a")
aplot::plot_list(gglist=list(alpha1,alpha2, alpha3, alpha4), tag_levels="a")
aplot::plot_list(gglist=list(beta1,beta2, beta3,beta4), tag_levels="a")
#aplot::plot_list(gglist=list(compo4,alpha4, beta4), tag_levels="a")


beta_compare <- ps_sites %>%
  mp_plot_dist(.distmethod = "bray", .group = Sites, group.test = TRUE,  textsize = 3) +
  scale_fill_manual(
    values = c("#9FBA95", "#f2b56f", "#7CAEF0", "#E6CECF", "#B696B6", "#80C1C4"),
    guide = guide_legend(keywidth = 0.6, keyheight = 0.6, label.theme = element_text(size = 10))
  ) +
  scale_color_manual(
    values = c("#9FBA95", "#f2b56f", "#7CAEF0", "#E6CECF", "#B696B6", "#80C1C4"),
    guide = guide_legend(keywidth = 0.6, keyheight = 0.6, label.theme = element_text(size = 10))
  ) +
  theme_classic() + 
  theme(panel.grid =element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 10,colour = "black"),
        legend.position = 'right')
  
beta_compare



p4 <- ps_sites %>% 
  mp_plot_ord(
    .ord = pcoa, 
    .group = Sites,
    #.starshape = Sites,
    .color = Sites, 
    show.adonis = T,
    #.size = Shannon, 
    #.alpha = Shannon,
    ellipse = T,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c("#9FBA95", "#E6CECF","#B696B6", "#80C1C4"), 
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
  ) +
  scale_color_manual(
    values=c("#9FBA95", "#E6CECF","#B696B6", "#80C1C4"),
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
  ) +
  scale_size_continuous(
    range=c(0.5, 3),
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
  ) 
p4
p5 <- ps_sites %>% 
  mp_plot_ord(
    .ord = pcoa, 
    .group = Places,
    .starshape = Sites,
    .color = Places, 
    show.adonis = T,
    #.size = Shannon, 
    #.alpha = Shannon,
    ellipse = F,
    show.legend = FALSE 
  ) +
  scale_fill_manual(
    values = c("#9FBA95", "#E6CECF","#B696B6", "#80C1C4"), 
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
  ) +
  scale_color_manual(
    values=c("#9FBA95", "#E6CECF","#B696B6", "#80C1C4"),
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
  ) +
  scale_size_continuous(
    range=c(0.5, 3),
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
  ) 
p5
aplot::plot_list(gglist=list(p4,p5), tag_levels="a")

#### Tree + diff
library(ggtree)
library(ggtreeExtra)
library(ggplot2)
library(MicrobiotaProcess)
library(tidytree)
library(ggstar)
library(forcats)
ps_sites %<>%
  mp_diff_analysis(
    .abundance = RelRareAbundanceBySample,
    .group = Sites,
    first.test.alpha = 0.01
  )

taxa.tree <- ps_sites %>% 
  mp_extract_tree(type="taxatree")
taxa.tree


f <- ps_sites %>%
  mp_plot_diff_cladogram(
    label.size = 2.5,
    hilight.alpha = .3,
    bg.tree.size = .5,
    bg.point.size = 2,
    bg.point.stroke = .25,
    taxa.class = "Family",
    as.tiplab = T
  ) +
  scale_fill_diff_cladogram( # set the color of different group.
    values = c("#9FBA95","#E6CECF","#B696B6", "#80C1C4")
  ) +
  scale_size_continuous(range = c(1, 4))
f

tree_all <- ggtree(
  taxa.tree,
  layout="radial",
  size = 0.3
) +
  geom_point(
    data = td_filter(!isTip),
    fill="white",
    size=1,
    shape=21
  ) +
  geom_hilight(
    data = td_filter(nodeClass == "Phylum"),
    mapping = aes(node = node, fill = label)
  ) +
  ggnewscale::new_scale("fill") +
  geom_fruit(
    data = td_unnest(RareAbundanceBySample),
    geom = geom_star,
    mapping = aes(
      x = fct_reorder(Sample, Sites, .fun=min),
      size = RelRareAbundanceBySample,
      fill = Sites,
      subset = RelRareAbundanceBySample > 0
    ),
    starshape = 13,
    starstroke = 0.25,
    offset = 0.04,
    pwidth = 0.8,
    grid.params = list(linetype=2)
  ) +
  scale_size_continuous(
    name="Relative Abundance (%)",
    range = c(.5, 3)
  ) +
  scale_fill_manual(values=c("#9FBA95", "#E6CECF","#B696B6", "#80C1C4"))+
  ggnewscale::new_scale("fill") +
  geom_fruit(
    geom = geom_col,
    mapping = aes(
      x = LDAmean,
      fill = Sign_Sites,
      subset = !is.na(LDAmean)
    ),
    orientation = "y",
    offset = 0.3,
    pwidth = 0.5,
    axis.params = list(axis = "x",
                       title = "Log10(LDA)",
                       title.height = 0.01,
                       title.size = 2,
                       text.size = 1.8,
                       vjust = 1),
    grid.params = list(linetype = 2)
  ) +
  ggnewscale::new_scale("size") +
  geom_point(
    data=td_filter(!is.na(Sign_Sites)),
    mapping = aes(size = -log10(fdr),
                  fill = Sign_Sites,
    ),
    shape = 21,
  ) +
  scale_size_continuous(range=c(1, 3)) +
  scale_fill_manual(values=c("#9FBA95", "#E6CECF","#B696B6", "#80C1C4"))+ theme(
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(0.3, "cm"),
    legend.spacing.y = unit(0.02, "cm"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 9),
  )+
  geom_tiplab(
    data = td_filter(!is.na(Sign_Sites) & nodeClass == "Genus"),
    aes(node = node, label = label),
    size = 5,
    offset = 0
  )

tree_all
p_map_plot = as.ggplot(p_map)
p2_plot = as.ggplot(p2)
p4_plot = as.ggplot(p4)

library(gridExtra)

layout <- rbind(c(1,2),
                c(1,2),
                c(3,2))              

grid.arrange(p_map_plot, p2_plot,p4_plot,layout_matrix = layout) 

p_map

p1

p2

tree_all


##### Input data of PVC####
setwd("S.frugi/1--Qiime")
# Entry 4 (Silva 138.1) on 1/28/2021 in R v3.6.3
S.ps <-qza_to_phyloseq(
  features="1--Data/table-with-phyla-no-metochondria-no-chloroplast.qza",
  tree="1--Data/rooted-tree.qza",
  taxonomy="1--Data/taxonomy.qza",
  metadata = "1--Data/metadata.tsv")

S.ps
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 727 taxa and 12 samples ]
#sample_data() Sample Data:       [ 12 samples by 3 sample variables ]
#tax_table()   Taxonomy Table:    [ 727 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 727 tips and 726 internal nodes ]
S.ps <- subset_samples(S.ps, Sites != "Excreted")
S.ps.dominate <- filter_taxa(S.ps, function(x) sum(x>2) > 1,TRUE)


T.ps <-qza_to_phyloseq(
  features="T.molitor/1--Data/otu_table.qza",
  tree="T.molitor/1--Data/rooted-tree.qza",
  taxonomy="T.molitor/1--Data/taxonomy.qza",
  metadata = "T.molitor/1--Data/group.txt")

T.ps
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 2389 taxa and 10 samples ]
#sample_data() Sample Data:       [ 10 samples by 1 sample variables ]
#tax_table()   Taxonomy Table:    [ 2389 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 2389 tips and 2384 internal nodes ]
#T.ps <- subset_samples(T.ps, Sites != "Excreted")
T.ps.dominate <- filter_taxa(T.ps, function(x) sum(x>2) > 2,TRUE)
T.ps.dominate

S_ps_M <- S.ps.dominate %>% as.MPSE() 
T_ps_M <- T.ps.dominate %>% as.MPSE() 

S_ps_M %<>%
  mp_cal_abundance( # for each samples
    .abundance = RareAbundance
  ) %>%
  mp_cal_abundance( # for each groups 
    .abundance=RareAbundance,
    .group=Intake
  )
T_ps_M %<>%
  mp_cal_abundance( # for each samples
    .abundance = RareAbundance
  ) %>%
  mp_cal_abundance( # for each groups 
    .abundance=RareAbundance,
    .group=Intake
  )

p1 <- S_ps_M %>%
  mp_plot_abundance(
    .abundance=RareAbundance,
    .group=Intake, 
    taxa.class = Genus, 
    topn = 20,
    relative = TRUE
  )
p1
p2 <- T_ps_M %>%
  mp_plot_abundance(
    .abundance=RareAbundance,
    .group=Intake, 
    taxa.class = Genus, 
    topn = 20,
    relative = TRUE
  )
p2


h1 <- S_ps_M %>%
  mp_plot_abundance(
    .abundance = RareAbundance,
    .group = Intake,
    taxa.class = Genus,
    relative = TRUE,
    topn = 20,
    geom = 'heatmap',
    features.dist = 'euclidean',
    features.hclust = 'average',
    sample.dist = 'bray',
    sample.hclust = 'average'
  )
h1


h2 <- T_ps_M %>%
  mp_plot_abundance(
    .abundance = RareAbundance,
    .group = Intake,
    taxa.class = Genus,
    relative = TRUE,
    topn = 20,
    geom = 'heatmap',
    features.dist = 'euclidean',
    features.hclust = 'average',
    sample.dist = 'bray',
    sample.hclust = 'average'
  )
h2
aplot::plot_list(gglist=list(h1, h2), tag_levels="A")

p3 <- S_ps_M %>%
  mp_plot_abundance(
    .abundance=RareAbundance, 
    .group=Intake,
    taxa.class = Genus,
    topn = 10,
    plot.group = TRUE
  )

p3
p4 <- T_ps_M %>%
  mp_plot_abundance(
    .abundance=RareAbundance, 
    .group=Intake,
    taxa.class = Genus,
    topn = 10,
    plot.group = TRUE
  )

p4

S_ps_M %<>% 
  mp_decostand(.abundance=Abundance)
S_ps_M %<>% mp_cal_dist(.abundance=hellinger, distmethod="bray")
T_ps_M %<>% 
  mp_decostand(.abundance=Abundance)
T_ps_M %<>% mp_cal_dist(.abundance=hellinger, distmethod="bray")

S_ps_M %<>% 
  mp_cal_pcoa(.abundance=hellinger, distmethod="bray")
T_ps_M %<>% 
  mp_cal_pcoa(.abundance=hellinger, distmethod="bray")


S_ps_M %<>%
  mp_adonis(.abundance=hellinger, .formula=~Intake, distmethod="bray", permutations=9999, action="add")
S_ps_M %>% mp_extract_internal_attr(name=adonis)

T_ps_M %<>%
  mp_adonis(.abundance=hellinger, .formula=~Intake, distmethod="bray", permutations=9999, action="add")
T_ps_M %>% mp_extract_internal_attr(name=adonis)


library(ggplot2)
### Alpha diversity
S_ps_M %<>% 
  mp_cal_alpha(.abundance=RareAbundance)
T_ps_M %<>% 
  mp_cal_alpha(.abundance=RareAbundance)

alpha1 <- S_ps_M %>% 
  mp_plot_alpha(
    .group=Intake, 
    test = "wilcox.test",
    .alpha=c(Chao1, ACE, Shannon, Simpson)
  ) +
  scale_fill_manual(values=c("#9FBA95", "#80C1C4"), guide="none") +
  scale_color_manual(values=c("#9FBA95",  "#80C1C4"), guide="none")
alpha2 <- T_ps_M %>% 
  mp_plot_alpha(
    .group=Intake, 
    test = "wilcox.test",
    .alpha=c(Chao1, ACE, Shannon, Simpson)
  ) +
  scale_fill_manual(values=c("#9FBA95","#80C1C4"), guide="none") +
  scale_color_manual(values=c("#9FBA95", "#80C1C4"), guide="none")
alpha1 + alpha2

beta1 <- S_ps_M %>% 
  mp_plot_ord(
    .ord = pcoa, 
    .group = Intake, 
    .color = Intake, 
    .size = Shannon, 
    show.adonis = T,
    #.alpha = Shannon,
    ellipse = TRUE,
    show.legend = FALSE # don't display the legend of stat_ellipse 
  ) +
  scale_fill_manual(
    values = c("#9FBA95",  "#80C1C4"), 
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
  ) +
  scale_color_manual(
    values=c("#9FBA95", "#80C1C4"),
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
  )+
  scale_size_continuous(
    range=c(2, 6),
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
  )

beta2 <- T_ps_M %>% 
  mp_plot_ord(
    .ord = pcoa, 
    .group = Intake, 
    .color = Intake, 
    .size = Shannon, 
    show.adonis = T,
    #.alpha = Shannon,
    ellipse = F,
    show.legend = FALSE # don't display the legend of stat_ellipse 
  ) +
  scale_fill_manual(
    values = c("#9FBA95",  "#80C1C4"), 
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
  ) +
  scale_color_manual(
    values=c("#9FBA95",  "#80C1C4"),
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
  )+
  scale_size_continuous(
    range=c(2, 6),
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
  )
beta1 / beta2

aplot::plot_list(gglist=list(h1, h2,beta1,beta2), tag_levels="a")


S_ps_M %<>%
  mp_diff_analysis(
    .abundance = RelRareAbundanceBySample,
    .group = Intake,
    first.test.alpha = 0.01
  )
T_ps_M %<>%
  mp_diff_analysis(
    .abundance = RelRareAbundanceBySample,
    .group = Intake,
    first.test.alpha = 0.01
  )

taxa.tree_S <- S_ps_M %>% 
  mp_extract_tree(type="taxatree",tip.level = "Genus")
taxa.tree_S

taxa.tree_T <- T_ps_M %>% 
  mp_extract_tree(type="taxatree",,tip.level = "Genus")
taxa.tree_T

phylum_colors <- c("Actinobacteriota" = "#9FBA95", "Bacteroidota" = "#E6CECF", 
                   "Firmicutes" = "#B696B6", "Proteobacteria" = "#80C1C4","p__un_d__Bacteria" = "grey")

tree_S <- ggtree(
  taxa.tree_S,
  layout="radial",
  size = 0.3
) +
  geom_point(
    data = td_filter(!isTip),
    fill="white",
    size=1,
    shape=21
  ) +
  geom_hilight(
    data = td_filter(nodeClass == "Phylum"),
    mapping = aes(node = node, fill = label)
  ) +
  scale_fill_manual(values = phylum_colors) + 
  ggnewscale::new_scale("fill") +
  geom_fruit(
    data = td_unnest(RareAbundanceBySample),
    geom = geom_star,
    mapping = aes(
      x = fct_reorder(Sample, Intake, .fun=min),
      size = RelRareAbundanceBySample,
      fill = Intake,
      subset = RelRareAbundanceBySample > 0
    ),
    starshape = 13,
    starstroke = 0.25,
    offset = 0.04,
    pwidth = 0.8,
    grid.params = list(linetype=2)
  ) +
  scale_size_continuous(
    name="Relative Abundance (%)",
    range = c(.5, 3)
  ) +
  scale_fill_manual(values=c("#9FBA95", "#80C1C4"))+
  ggnewscale::new_scale("fill")  + geom_tiplab(size=3, offset=7.2)

tree_S

phylum1_colors <- c("Actinobacteriota" = "#9FBA95", "Bacteroidota" = "#E6CECF", "Firmicutes" = "#B696B6"
                    , "Proteobacteria" = "#80C1C4","Chloroflexi" = "#fae69e","Dadabacteria" = "#f2b56f",
                    "Fusobacteria" = "#BB9393","Patescibacteria" = "#808080","Spirochaetes" = "#7CAEF0",
                    "Tenericutes" = "#cee9dc","Verrucomicrobia" = "#ffadbb")

taxa.tree_T@phylo$tip.label <- sub("^g__", "", taxa.tree_T@phylo$tip.label)
taxa.tree_T@phylo$node.label<- sub("^p__", "", taxa.tree_T@phylo$node.label)

tree_T <- ggtree(
  taxa.tree_T,
  layout="radial",
  size = 0.3
) +
  geom_point(
    data = td_filter(!isTip),
    fill="white",
    size=1,
    shape=21
  ) +
  geom_hilight(
    data = td_filter(nodeClass == "Phylum"),
    mapping = aes(node = node, fill = label)
  ) +
  scale_fill_manual(values = phylum1_colors) + 
  ggnewscale::new_scale("fill") +
  geom_fruit(
    data = td_unnest(RareAbundanceBySample),
    geom = geom_star,
    mapping = aes(
      x = fct_reorder(Sample, Intake, .fun=min),
      size = RelRareAbundanceBySample,
      fill = Intake,
      subset = RelRareAbundanceBySample > 0
    ),
    starshape = 13,
    starstroke = 0.25,
    offset = 0.04,
    pwidth = 0.8,
    grid.params = list(linetype=2)
  ) +
  scale_size_continuous(
    name="Relative Abundance (%)",
    range = c(.5, 3)
  ) +
  scale_fill_manual(values=c("#9FBA95", "#80C1C4"))+
  ggnewscale::new_scale("fill") + geom_tiplab(size=3, offset=7.2)

tree_T

p1 + p2
p3 + p4
alpha1 + alpha2 + beta1 +beta2
tree_S + tree_T

aplot::plot_list(gglist=list(beta1, beta2,tree_S,tree_T), tag_levels="a")


#### Differential Enzyme

library(pipeR)
library(pheatmap)
#setwd("Picrust_gene")
#setwd("Picrust2/picrust2_result_pipeline2/EC_metagenome_out/Picrust_gene")
otu <- read.delim("Ec_description_PVC1.csv", header=T, row.names=1,sep = ",")
otu$EF <- NULL
otu$PEF <- NULL

otu=as.data.frame(((otu)/rowSums(otu)*100))

first <- colorRampPalette(c("#9FBA95", "white"))(50)
second <- colorRampPalette(c("white", "#80C1C4"))(50)
palette <- c(first, second)
p2 <- pheatmap(otu,treeheight_col=5,color=palette,
               cellwidth=12,  cellheight=10,
               #annotation_col=annot_data,
               cluster_cols =F,
               #display_numbers = T,
               cluster_row=T,show_rownames=T, width=4, height=4,main='Plastic degradation gene relative richness')
wilcox.test(otu)

plot1 <- as.ggplot(p1)
plot2 <- as.ggplot(p2)

plot1 + plot2


#### Differential Pathway T. molitor
library(DESeq2)
setwd("Picrust2/picrust2_result_pipeline2/KO_metagenome_out/Deseq2")
mycounts <- read.csv("KO_IN.csv")
head(mycounts)
rownames(mycounts)<-mycounts[,1]
mycounts<-mycounts[,-1]
head(mycounts)
condition <- factor(c(rep("Normal",4),rep("PVC",6)), levels = c("Normal","PVC"))
condition
colData <- data.frame(row.names=colnames(mycounts), condition)
colData
dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
dds_norm <- DESeq(dds)
dds_norm
normalized_counts <- counts(dds_norm, normalized=TRUE)
head(normalized_counts)

normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
write.table(normalized_counts, file="dds_normalized_counts.xls",
            quote=F, sep="\t", row.names=T, col.names=T)
system(paste("sed -i '1 s/^/ID\t/'", "dds_normalized_counts.xls"))
rld <- rlog(dds_norm, blind=FALSE)
rlogMat <- assay(rld)
rlogMat <- rlogMat[order(normalized_counts_mad, decreasing=T), ]
write.table(rlogMat, file="dds_normalized_counts_rlog.xls",
            quote=F, sep="\t", row.names=T, col.names=T)
system(paste("sed -i '1 s/^/ID\t/'", "dds_normalized_counts_rlog.xls"))

res = results(dds_norm, contrast=c("condition", "PVC", "Normal"))
res = res[order(res$pvalue),]
head(res)
summary(res)
write.csv(res,file="All_results_OUT.csv")
table(res$padj<0.01)
diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
dim(diff_gene_deseq2)
head(diff_gene_deseq2)
write.csv(diff_gene_deseq2,file= "DEG_p<0.05_log>1.csv") 

library(ggplot2)
library(BiocGenerics)
library(EnhancedVolcano)
library(airway)
p2 <- EnhancedVolcano(res,
                      lab = rownames(res),
                      x = 'log2FoldChange',
                      y = 'padj',
                      xlim = c(-12, 10),
                      title = 'Corn versus PVC',
                      pCutoff = 0.01,
                      FCcutoff = 2,
                      col=c('grey', '#9FBA95', 'black', '#80C1C4'),
                      colAlpha = 0.5,
                      #legend=c('NS','log2FoldChange','P value',
                      #         'P value & log2FoldChange'),
                      legendPosition = 'right',
                      legendLabSize = 14,
                      legendIconSize = 5.0,
                      selectLab = character(0)
)
p2
#### Differential Pathway S. frugiperda

library(DESeq2)
setwd("5-Picrust/picrust2_result_pipeline2/KO_metagenome_out")
mycounts <- read.table('pred_metagenome_unstrat.tsv', header=T, sep='\t')
head(mycounts)
rownames(mycounts)<-mycounts[,1]
mycounts<-mycounts[,-1]
head(mycounts)
mycounts <- round(mycounts)
mycounts <- mycounts[,c("IN.IF1","IN.IF2","IN.IF3","IN.PIF1","IN.PIF2","IN.PIF3")]
condition <- factor(c(rep("Normal",3),rep("PVC",3)), levels = c("Normal","PVC"))
condition
colData <- data.frame(row.names=colnames(mycounts), condition)
colData
dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
dds_norm <- DESeq(dds)
dds_norm
normalized_counts <- counts(dds_norm, normalized=TRUE)
head(normalized_counts)

normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
write.table(normalized_counts, file="dds_normalized_counts.xls",
            quote=F, sep="\t", row.names=T, col.names=T)
system(paste("sed -i '1 s/^/ID\t/'", "dds_normalized_counts.xls"))
rld <- rlog(dds_norm, blind=FALSE)
rlogMat <- assay(rld)
rlogMat <- rlogMat[order(normalized_counts_mad, decreasing=T), ]
write.table(rlogMat, file="dds_normalized_counts_rlog.xls",
            quote=F, sep="\t", row.names=T, col.names=T)
system(paste("sed -i '1 s/^/ID\t/'", "dds_normalized_counts_rlog.xls"))

res = results(dds_norm, contrast=c("condition", "PVC", "Normal"))
res = res[order(res$pvalue),]
head(res)
summary(res)
write.csv(res,file="All_results_OUT.csv")

table(res$padj<0.01)
diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
dim(diff_gene_deseq2)
head(diff_gene_deseq2)
write.csv(diff_gene_deseq2,file= "DEG_p<0.01_log>2.csv") 

library(ggplot2)
library(BiocGenerics)

#BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
library(airway)
p1 <- EnhancedVolcano(res,
                      lab = rownames(res),
                      x = 'log2FoldChange',
                      y = 'padj',
                      xlim = c(-12, 10),
                      title = 'Corn versus PVC',
                      pCutoff = 0.01,
                      FCcutoff = 2,
                      col=c('grey', '#9FBA95', 'black', '#80C1C4'),
                      colAlpha = 0.5,
                      #legend=c('NS','log2FoldChange','P value',
                      #         'P value & log2FoldChange'),
                      legendPosition = 'right',
                      legendLabSize = 14,
                      legendIconSize = 5.0,
                      selectLab = character(0)
)
p1



plot1 <- as.ggplot(p1)
plot2 <- as.ggplot(p2)
plot1 + plot2


#### China Map#### 
library(ggplot2)
library(sf)

china_map <- st_read("中华人民共和国.shp")
china_map$name

provinces <- c("云南省", "海南省", "河南省", "湖北省")

ggplot(data = china_map) +
  geom_sf(aes(fill = name %in% provinces)) +
  scale_fill_manual(values = c("TRUE" = "orange", "FALSE" = "grey")) +
  labs(title = "Map of China Highlighting Yunnan, Hainan, Henan, and Hubei") +
  theme_minimal()

province_colors <- c("云南省" = "#80C1C4", "海南省" = "#B696B6", "河南省" = "#E6CECF", "湖北省" = "#9FBA95")

p_map <- ggplot(data = china_map) +
  geom_sf(color = "white", size = 0.2, aes(fill = ifelse(name %in% provinces, name, "Other"))) +
  scale_fill_manual(
    values = c(province_colors, "Other" = "#BDC3C7"),
    labels = c("Yunnan", "Hainan", "Henan", "Hubei", "Other Provinces")
  ) +
  labs(
    title = "Map of China Highlighting Yunnan, Hainan, Henan, and Hubei",
    subtitle = "Each highlighted province is shown in a unique color",
    caption = "Source: GADM"
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12),
    plot.caption = element_text(size = 10, color = "grey50"),
    legend.position = "none"#,  # Put the legend to the right
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    #panel.background = element_rect(fill = "#F3F4F6", color = NA),
    #plot.background = element_rect(fill = "#F3F4F6", color = NA)
  )
p_map


