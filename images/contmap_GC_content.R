library(ggtree)
library(cowplot)
library(phytools)
library(tidyverse)
library(viridis)
library(patchwork)
library(gridExtra)
setwd("C:\\Users\\johns\\OneDrive - University of Calgary\\data\\projects\\esc2025_comp_gen\\response")
AAtree<-read.tree('trimmed.ms.lg.treefile')
NTtree<-read.tree('nt3.treefile')
out<-c(AAtree$tip.label[grepl('Chao', AAtree$tip.label)], AAtree$tip.label[grepl('CLUM', AAtree$tip.label)],
       AAtree$tip.label[grepl('Pvander', AAtree$tip.label)], AAtree$tip.label[grepl('Moch', AAtree$tip.label)],
       AAtree$tip.label[grepl('Culic', AAtree$tip.label)])

out

AAroot<-root(phy = AAtree, resolve.root = TRUE, outgroup = out)
NTroot<-root(phy = NTtree, resolve.root = TRUE, outgroup = out)


dat<-read.table('gc.tsv')
dat<-read.table('gc.txt', header = T, sep = '\t')
dat.plot<-data.frame(tax = dat$Sequence_ID, gc = dat$GC, pos = dat$pos) %>%
  # mutate(pos2 = case_when(pos == 1 | pos == 2 ~ 'pos12', .default = 'pos3')) %>%
  mutate(position = case_when(pos == 1 ~ 'pos1',
                              pos == 2 ~ 'pos2',
                              pos == 3 ~ 'pos3')) %>%
  group_by(tax, position) %>%
  summarise(gc = mean(gc)) %>%
  drop_na %>%
  pivot_wider(names_from = position, values_from = gc)

####### Boxplots
is_tip <- AAroot$edge[,2] <= length(AAroot$tip.label)
ordered_tips <- AAroot$edge[is_tip, 2]

level_order <- AAroot$tip.label[ordered_tips]
p2<-ggplot(dat,aes(x=GC,y=factor(Sequence_ID,level=level_order),fill=as.factor(pos)))+
  geom_boxplot(show.legend=FALSE,outliers=F,median.color = "white")+
  scale_fill_viridis(discrete=TRUE)+
  theme_bw()+ylab("Taxa")+
  xlab("GC")+facet_wrap(vars(pos),nrow=1)
p2
#theme(axis.text.y = element_blank())+
grid.arrange(p,p2,nrow=1)

mansoniini2_levels <- c("FMEL09", "FMEL10", "YF11", "Mnsept_Au50")

outgroup_levels <- c("Chaoborus_flavidulus_GGBK01_1.1",
                     "CLUMA_GCA_900005825_CLUMA1",
                     "Pvanderplanki_20121018",
                     "Mochlonyx_cinctipes_GCA_001014845_1",
                     "Culicoides_sonorensis_GCA_900002565_1")

# Add new grouping column
df <- dat %>%
  mutate(Group = case_when(
    Sequence_ID %in% mansoniini2_levels ~ "Mansoniini2",
    Sequence_ID %in% outgroup_levels ~ "Outgroups",
    TRUE ~ "Ingroup"   # everything else
  )) %>%
  mutate(Group = factor(Group, levels = c("Ingroup", "Mansoniini2", "Outgroups")))

ggplot(df,aes(x=GC,y=factor(Group),fill=as.factor(Position)))+
  geom_boxplot(show.legend=FALSE,outliers=F,median.color = "white")+
  scale_fill_viridis(discrete=TRUE)+
  theme_bw()+ylab("Taxa")+
  xlab("GC")+facet_wrap(vars(Position),nrow=1)

ggplot(df,aes(x=GC,y=factor(Group),fill=as.factor(Position)))+
  geom_boxplot(show.legend=FALSE,outliers=F,median.color = "white")+
  scale_fill_viridis(discrete=TRUE)+
  theme_bw()+ylab("Taxa")+
  xlab("GC")

gc1<-setNames(dat.plot$pos1, dat.plot$tax)
gc2<-setNames(dat.plot$pos2, dat.plot$tax)
gc3<-setNames(dat.plot$pos3, dat.plot$tax)


############with the amino acid tree
fitAA1<-fastAnc(AAroot, gc1, vars = TRUE, CI = TRUE)
fitAA2<-fastAnc(AAroot, gc2, vars = TRUE, CI = TRUE)
fitAA3<-fastAnc(AAroot, gc3, vars = TRUE, CI = TRUE)


AAtd1<-data.frame(node = nodeid(AAroot, names(gc1)),
                  trait = gc1)
AAtd2<-data.frame(node = nodeid(AAroot, names(gc2)),
                  trait = gc2)
AAtd3<-data.frame(node = nodeid(AAroot, names(gc3)),
                  trait = gc3)
AAnd1<-data.frame(node = names(fitAA1$ace), trait = fitAA1$ace)
AAnd2<-data.frame(node = names(fitAA2$ace), trait = fitAA2$ace)
AAnd3<-data.frame(node = names(fitAA3$ace), trait = fitAA3$ace)

AA1<-rbind(AAtd1,AAnd1)
AA2<-rbind(AAtd2,AAnd2)
AA3<-rbind(AAtd3,AAnd3)

AA1$node<-as.numeric(AA1$node)
AA2$node<-as.numeric(AA2$node)
AA3$node<-as.numeric(AA3$node)

AA1root2<-full_join(AAroot, AA1, by = 'node')
AA2root2<-full_join(AAroot, AA2, by = 'node')
AA3root2<-full_join(AAroot, AA3, by = 'node')
p <- ggtree(AAroot,ladderize=FALSE)
facet_plot(p, panel = "Sequence Distance", data = dat.plot_long, geom_boxplot, 
           mapping = aes(x=value,y=as.factor(pos),fill = as.factor(pos)), alpha = .6)


(aa1.plot<-ggtree(AA1root2,
                  aes(color = trait),
                  continuous = 'color', linewidth  = 0.6) +
    # xlim(c(0, 1.5)) +
    # scale_color_viridis('GC-content (%)') +
    scale_color_viridis('GC-content (%)',limits = c(0.1, 0.8)) +
    ggtitle('Nucleotide position 1') +
    # geom_tiplab(color = 'black', size = 3) +
    geom_treescale(x = 0.05, y= 32) )

(aa2.plot<-ggtree(AA2root2,
                  aes(color = trait),
                  continuous = 'color', linewidth  = 0.6) +
    # xlim(c(0, 1.5)) +
    # scale_color_viridis('GC-content (%)') +
    scale_color_viridis('GC-content (%)',limits = c(0.1, 0.8)) +
    ggtitle('Nucleotide position 2') +
    # geom_tiplab(color = 'black', size = 3) +
    geom_treescale(x = 0.05, y= 32) )
#oob=scales::squish,

(aa3.plot<-ggtree(AA3root2,
                  aes(color = trait),
                  continuous = 'color', linewidth  = 0.6) +
    # xlim(c(0, 1.5)) +
    # scale_color_viridis('GC-content (%)') +
    scale_color_viridis('GC-content (%)',limits = c(0.1, 0.8)) +
    ggtitle('Nucleotide position 3') +
    # geom_tiplab(color = 'black', size = 3) +
    geom_treescale(x = 0.05, y= 32) )


(panel1<-aa1.plot + aa2.plot + aa3.plot + plot_layout(guides = 'collect'))

############
############
############
############
############
############
############
############

########### with the NT3 tree
fitNT1<-fastAnc(NTroot, gc1, vars = TRUE, CI = TRUE)
fitNT2<-fastAnc(NTroot, gc2, vars = TRUE, CI = TRUE)
fitNT3<-fastAnc(NTroot, gc3, vars = TRUE, CI = TRUE)


NTtd1<-data.frame(node = nodeid(NTroot, names(gc1)),
                  trait = gc1)
NTtd2<-data.frame(node = nodeid(NTroot, names(gc2)),
                  trait = gc2)
NTtd3<-data.frame(node = nodeid(NTroot, names(gc3)),
                  trait = gc3)
NTnd1<-data.frame(node = names(fitNT1$ace), trait = fitNT1$ace)
NTnd2<-data.frame(node = names(fitNT2$ace), trait = fitNT2$ace)
NTnd3<-data.frame(node = names(fitNT3$ace), trait = fitNT3$ace)

NT1<-rbind(NTtd1,NTnd1)
NT2<-rbind(NTtd2,NTnd2)
NT3<-rbind(NTtd3,NTnd3)

NT1$node<-as.numeric(NT1$node)
NT2$node<-as.numeric(NT2$node)
NT3$node<-as.numeric(NT3$node)

NT1root2<-full_join(NTroot, NT1, by = 'node')
NT2root2<-full_join(NTroot, NT2, by = 'node')
NT3root2<-full_join(NTroot, NT3, by = 'node')


(nt1.plot<-ggtree(NT1root2,
                  aes(color = trait),
                  continuous = 'color', size  = 0.6) +
    # xlim(c(0, 1.5)) +
    # scale_color_viridis('GC-content (%)') +
    scale_color_viridis('GC-content (%)', limits = c(0.1, 0.8)) +
    ggtitle('Nucleotide position 1') +
    # geom_tiplab(color = 'black', size = 3) +
    geom_treescale(x = 0.05, y= 32) )

(nt2.plot<-ggtree(NT2root2,
                  aes(color = trait),
                  continuous = 'color', size  = 0.6) +
    # xlim(c(0, 1.5)) +
    # scale_color_viridis('GC-content (%)') +
    scale_color_viridis('GC-content (%)', limits = c(0.1, 0.8)) +
    ggtitle('Nucleotide position 2') +
    # geom_tiplab(color = 'black', size = 3) +
    geom_treescale(x = 0.05, y= 32) )

(nt3.plot<-ggtree(NT3root2,
                  aes(color = trait),
                  continuous = 'color', size  = 0.6) +
    # xlim(c(0, 1.5)) +
    # scale_color_viridis('GC-content (%)') +
    scale_color_viridis('GC-content (%)', limits = c(0.1, 0.8)) +
    ggtitle('Nucleotide position 3') +
    #geom_tiplab(color = 'black', size = 3) +
    geom_treescale(x = 0.05, y= 32) )

panel2<-nt1.plot + nt2.plot + nt3.plot + plot_layout(guides = 'collect')


############
############
############
############
############save plots
ggsave(plot = panel1, filename = 'AA_contmap_panel2.pdf', width = 8, height = 8)
ggsave(plot = panel2, filename = 'NT_contmap_panel2.pdf', width = 8, height = 8)



library(ape)
library(ggtree)
library(ggplot2)

# Random tree with 10 tips
tree <- rtree(10)

# Multiple continuous variables for each tip
tip_data <- data.frame(
  label = tree$tip.label,
  trait1 = runif(10, 0, 1),
  trait2 = runif(10, 5, 10),
  trait3 = rnorm(10, 50, 5)
)

head(tip_data)
head(dat.plot)
library(tidyr)

dat.plot_long <- dat.plot %>%
  pivot_longer(cols = starts_with("pos"),
               names_to = "pos",
               values_to = "value")
head(dat.plot_long)

library(ape)
library(ggtree)
library(ggplot2)

# Random tree with 10 tips
tree <- rtree(10)

# Continuous variables for each tip
tip_data <- data.frame(
  label = tree$tip.label,
  trait1 = runif(10, 0, 1),
  trait2 = runif(10, 5, 10),
  trait3 = rnorm(10, 50, 5)
)
library(tidyr)
tip_data_long <- tip_data %>%
  pivot_longer(cols = starts_with("trait"),
               names_to = "trait",
               values_to = "value") %>%
  mutate(trait = factor(trait),
         trait_idx = as.integer(trait))  # numeric x required by facet_plot

# Base tree
p <- ggtree(tree)

# Add heatmap panel with facet_plot (use numeric x)
g <- facet_plot(p, panel = "heatmap", data = tip_data_long,
                geom = geom_tile,
                mapping = aes(x = trait_idx, y = label, fill = value),
                width = 0.9, height = 0.9)
g
g +
  scale_x_continuous(breaks = levels(tip_data_long$trait) |> as.integer(),
                     labels = levels(tip_data_long$trait)) +
  scale_fill_viridis_c(option = "magma") +
  theme_minimal() +
  theme(
    panel.spacing.x = unit(0.5, "lines"),
    strip.text.x = element_text(size = 10),
    axis.title = element_blank(),
    axis.text.y = element_blank()
  )
