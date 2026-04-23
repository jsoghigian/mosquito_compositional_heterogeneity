
#Introduction
by John Soghigian, Gen Morinaga, and Huiqing Yeo
Recently, colleagues and I were interested in investigating patterns in multiple sequence alignments that resulted in disagreement between amino acids and nucleotide position three, which resulted in topological differences based primarily on the position of the root. For the particular organismal group of interest, mosquitoes, the two topologies were one in which the two mosquito subfamilies were monophyletic and the root was on the branch between them, or another in which the root was within one of the subfamilies (The Culicinae), resulting in the subfamily being non-monophyletic, a pattern described in Pierce et al. 2025. Based on a range of past evidence, we had strong a priori assumptions about which of these was likely to be the case. We felt that it was likely that compositional heterogeneity, a well-known problem in maximum likelihood phylogenetic inference, could create such a conflict. We decided to delve, in detail, into the causes - and their consequences - of long branches and compositional heterogeineity.

# Datasets
Here, we analyze four datasets, from previously published manuscripts, with a minor exception (Dataset B uses different outgroups, which are publicly available). We could unfortunately not analze the datasets from Pierce et al., as their alignments were not released. None the less, Dataset A mimics their dataset quite well, including an over-sampling of Anophelines, relatively low sampling of tribes, and choice of GC-poor outgroups.
-Dataset A: A dataset previously published in Soghigian et al. 2023, of 5667 genes from 54 species of mosquitoes and outgroup Diptera, based on genomes and transcriptomes. This is the main dataset we will work with.
-Dataset B: A dataset previously published in Soghigian et al. 2023, of 5667 genes from 54 species of mosquitoes but with new, publicly available outgroups from genomes and transcriptomes. This dataset is just used briefly to demonstrate how choice of outgroup, in particular as it relates to GC content, changes root position substantially.
-Dataset C: A dataset also published in Soghigian et al. 2023, of 709 genes from 260 mosquitoes and outgroup Diptera, based on anchor hybrid enrichment, genomes, and transcriptomes. We do show some results from this dataset, but most of our analyses focus on Dataset A.
-Dataset D: A dataset of UCEs derived from a genome dataset previously published in Morinaga et al. 2025, of 11 mosquito genomes and 3 outgroup Diptera genomes. This dataset was used primarily to characterize UCEs.

# Phylogenetic Relationships in Mosquitoes, from Amino Acids and NT3 Alone
Datasets A and C are protein-coding, which enables us to analyze amino acids as well as any codon position jointly or separately. Previous analyses on these datasets indicate a phylogeny similar to below, in which the two subfamilies were monophyletic. These analyses were based on amino acids and codon position 2 or 1 and 2 together (hereafter NT1 or NT1/2), and resolved the same phylogeny. Note unsampled tribes indicated by dashed lines.
<Figure 1 here with Caption>
Codon position 3 (hereafter NT3) had been excluded due to evidence of significant saturation (See Soghigian et al. 2023 supplemental). However, a recent paper argued that NT3 was a more reliable phylogenetic marker than amino acids or NT1/2. Although virtually all evidence outside of codon position 3 and UCEs support two monophyletic subfamilies, a recent publication proposed a topology based predominantly on signal from NT3. We can mimic a similar signal in our datasets A and C by analyzing position 3:
```
iqtree3 -s nt3.fasta -spp partitions.pos3.nex -pre nt3.MFP -m MFP -B 1000 -nt 16
```
If you are wanting to do this on your own data, you could either use a program like AMAS to extract position 3 given a relevant partition file, or you could also just specify the partition file to IQTree in a way that it only analyzes position three. For instance:
```
DNA, OG10_1 = 1-13491\3
DNA, OG10_2 = 2-13491\3
DNA, OG10_3 = 3-13491\3
```
Will cause IQ-Tree to analyze each codon position separately. If you were to delete those first two lines, you would analze only NT3.
Whether using ModelFinder (MFP) or, say, GTR+F+I+R, the topology reconstructed from NT3 alone results in a non-monophyletic Culicinae, as below:
<Figure with NT3 topologies)
Now both topologies do have some problematic aspects, as tribes are not monophyletic. Note that this same phenomenon was observed by Pierce et al. in their analyses, where the Mansoniniini was not monophyletic as Mansonia and Coquelltidia were in separate parts of the phylogeny. This is not a problem we see with amino acids or other nucleotide positions.

But why does NT3 show this topology? That's the question

## Assessing Support for Conflicting Topologies
There are a variety of ways to assess support for topologies that conflict. One option is to assess the difference in log-likelihood between two topologies given an alignment. 
<Images>

That likelihood can be assesed either by site, or by partition. Estimating the likelihood can be done as below:
```
iqtree3 -s gen_tran.trim.dna.phy -spp pos123.lnl.best_model.nex -z two_topos.trees -wsl -wpl -n 0 -nt 16 --prefix pos123.lnl
```
With this command, IQ-Tree will output the log likelihood by site (-wsl) and by partition (-wpl). As we are primarily interested in conflict at the third codon position, we specify a partition scheme that separates each gene by codon position and we then can analyze the resulting partition log-likeilhood file. IQ-Tree outputs a .partlh file following the above command, with one row for each tree supplied with -z. Each partition is separated by a space, and the order of partitions is dictated by the order of partitions.

As we run these sort of analyses often, we've created a simple tool to merge together the .iqtree and lh files:
```
lnleval.py file.iqtree lh.partlh
```
The resulting tsv can be imported into something like R and processed and summarized. We've uploaded a file that includes these values, ntaastats.tsv. If one summarzies the support across different partition types from Dataset A, we find the following for this dataset:

This aligns with expectations given the topologies we recover from amino acids and other positions. But again - why?

## Inspecting Alignments
There has been significant prior work describing what characteristics are favorable or unfavorable in multisequence alignments, such as calculations of saturation, evolutionary rate, etc. To evaluate NT3 in detail, we utilized the functions in the excellent software package PhyKit. Another option is GeneSortR, which can produce many of the same statistics en masse, but requires a stricter input format.

We decided to evalute each partition type (and amino acids, at least for many statistics) for saturation  (lower values are usually desierable), mean patristic distance (lower values are usually desierable), evolutionary rate (lower values are usually desierable), treeness/relative compositional variance (higher values usually desierable), long-branch score (higher values are often associated with long branch artifacts), and GC-content (lower values are usually desirable, though this is taxon-specific). For the last two, we considered not only each nucleotide codon position in a partition, but we also considered per-taxon measures. When using PhyKit, each alignment and/or tree must be specified individually; we collated results together using simple grep commands in bash. Consult PhyKit documents for specific commands relative to each of the measures we referred to above.
```
cd aln/
for file in $(ls *.fasta);
do
og=$(basename $file .fasta)
phykit gc_content $file -v > gctmp.txt
paste gctmp.txt <(for line in $(awk -F'\t' '{print $2}' gctmp.txt);do echo ${og};done) >> ../gc.txt
rm gctmp.txt
done
```

Where commands required both the trees and the alignment, we ran loops such as this:
```
cd trees/
for file in $(ls *.treefile);
do
og=$(basename $file .treefile)
sat_val=$(phykit saturation -t $file -a ../aln/$og.fasta)
echo $og $sat_val >> ../sat_og.txt
unset sat_val
unset og
done
```
After calculating values across an entire alignment, we merged together various parameters into a single TSV. We used Excel for this, but it can also be accomplished (perhaps more simply) in R. We then merged the log likelihood values from analyzing each partition via IQ-Tree3 and the result was the attached datafile ntaastats.tsv.

Plots for Saturation, Average Patristic Distance, Evolutionary Rate, and Treeness/RCV

Unsururpsingly based on expectations of third codon position, our calculations indicate that the third codon position has evidence of significant saturation (consistent with Soghigian et al. 2023), higher average patristic distances, a faster evolutionary rate, and lower treeness/RCV, all of which indicate this partition likely has lower phylogenetic value than other partitions.

We decided to explore long branch score and GC content on both a per partition and a per taxon level. This was for two reasons: 1) the topology for NT3 suggested the possibility of long-branch attraction for some species, and 2) Pierce et al. 2025 noted some differences in GC content in the third codon position which appeared to suggest compositional heterogeneity at this position. High levels of compositional heterogeneity in phylogenetics is a major source of error in maximum likelihood inference, and one of the assumptions of ML inference is compositional homogeneity of sequences in the alignment. 

Two Plots for Long Branch Score and GC Content Per Taxa

Both of these measures confirmed our suspicions: Some taxa had much longer branches for trees built from the third codon position than they did for the first two, and GC content was highly heterogeneous at the third codon position. Thus, it appeared that saturation at NT3, as well as codon usage variability between species, was contributing to long branch artifacts and compositional heterogeineity - both of which could result in significant errors in phylogenetic inference. 

We can plot GC content over the topologies to observe this heterogeneity in another way, using the R package phytools and its funtion contmap (see contmaps2.R for code). 


Here, we can clearly see that NT1 and NT2 have limited compositional heterogeneity between branches, while NT3 has considerable variability. The colored branches across the phylogeny are an ancestral state reconstruction; internal branch color values should be treated with caution.


## Ameloirating Compositional Heterogeneity
Having identified that NT3 had substantial evidence of saturation and compositional heterogeineity, we turned to the literature to address this. There are a range of methods that have been developed to address these issues, such as mixture models (which have been shown in amino acid datasets to resolve long branch attraction issues to some extent), and recoding. Mixture models have been implemented in various software suites, such as PhyloBayes and IQ-Tree. Recoding involves a strategy in which nucleotides or amino-acids are recoded to a reduced set of character states meant to reflect a degree of exchangeability based around either substitution models or known phenomenon. Mixture models are particularly computational intensive, and we did not observe a difference in topology when using them in the past, so we focused on how to handle the NT3 data.

For NT3, a natural choice is to use the RY recoding scheme, wherein purines (A and G) are recoded to R and pyrimidines (T and C) are recoded to Y. RY recoding has previously been used in diverse branches of the tree of life to ameloirate compositional heterogeineity, especially at the third codon position. There are more complex nucleotide recoding schemes (e.g., see Noah et al. 2020 and 'principled recoding' around degenerate codons), but simple RY seemed a good place to start.

In addition, to evaluate if the amino acid topologies were a result of long branch attraction, we also recoded the amino acid dataset following a Dayhoff 6, Dayhoff 15, and Dayhoff 18 strategy.

All recoding was conducted using PhyKit, e.g:
```
phykit alignment_recoding nt3.fasta RY-recoding.txt
```

We recoded both Dataset A and Dataset C third codon position and ran IQ-Tree3 on the resulting recoded matrix. We ran analyses where we treated these alignments either as binary, or as nucleotide data, with ModelFinder Plus. We chose to evaluate the binary option as there have been suggestions to analyze recoded data using non-specific substitution models originally created for morphology (see [here](https://iqtree.github.io/doc/Substitution-Models#binary-and-morphological-models)). Regardless, the choice of treating data as binary or note did not make a difference in the topology.

<TREE FILES FOR NT3 RECODED)

Recoding the 3rd codon position in both Dataset A and B resulted in strong support for both subfamilies as monophyletic. All branch support values for all nodes at or above the tribe level were 99 or 100 (available in the logs_trees directory).

In fact, recoding all 3 codon positions results in strong overall support for the monophyletic Culicinae:

<BAR PLOTS GEN MADE FOR RECODED>

What about recoding amino acids? To expedite tree inference, ModelFinderPlus was not used, but a GTRX+FO model was used that is a general time reversible model with unequal rates for each state, with state frequences estimated from the data. Regardless of recoding scheme, recoded the amino acids reflected essentially the same higher-level topology as the unrecoded data. As such, we only present the Dayhoff 15 recoded data, below:

DAYHOFF

This provides strong evidence that compositional heterogeineity was likely to blame for an incorrect root placement in NT3 data, likely owed to long branch attraction.

## Simulating Alignments
One might wonder if alignment recoding itself could introduce errors. There are cases where alignment recoding does reduce the accuracy of phylogenetic inference (see [Foster et al. 2023](https://doi.org/10.1093/sysbio/syac042) for some details on recoding), and in general, it is best to interpret trees from recoded data with some degree of caution. In our case, however, recoding NT3 resolves the discrepancy between this position and other positions/amino acids - precisely as expected if compositional heterogeineity was at fault for the topological discordance.

As a further sanity check, we can also simulate NT3 alignments under the discordant NT3 topology, based on the alignment and overall GC content of NT3 (but without compositional heterogeineity), and evaluate if alignment recoding influences the topology that results. If it did, this might suggest that recoding itself may be creating a phylogenetic artifact, rather than resolving one.

Simulating alignments is now extremely easy with IQ-Tree.  See below:
```
iqtree3 -s nt3.fasta -p partitions.pos3.nex --alisim mimicked_nt3_gtr -m GTR+F+I+R --seqtype DNA -t nt3.treefile -nt 16 --out-format fasta --num-alignments 100
```

This command simualted 100 alignments based on the alignment and partition scheme of the nt3 fasta file, under the model GTR+F+I+R, with the same base frequencies as the input alignment but simulated alignments were compositionally homogenous. Unsurprisingly, all of these alignments, when recoded, recover the same topology as unrecoded data. We could evaluate each alignment's log-likelihood support for either the original two subfamilies or the NT3 topology, or we could compare distances between trees, or any other such metric. When topologies are relatively similar, one might also consider something like DensiTree. In our case, the topologies didn't change due to recoding at all.

So... it looks like the NT3 data from this dataset do indeed support two monophyletic subfamilies!

## Choice of Outgroups

As previously mentioned, alignent recoding could reduce the accuracy of phylogenetic inference, and there are some strong feelings against it. Although we feel that its usage is appropriate to evaluate if compostional heterogeneity is responsible for topological artifacts, another option might be to address the problem this dataset has head on: Namely, that outgroups have significantly lower GC content at third position than most ingroups. If the GC content of outgroups was causing such significant compositional heterogeneity that it interfered with phylogenetic inference, then swapping the lowest GC content outgroups for higher GC content outgroups might resolve this.

The lowest GC content outgroups are Clunio marinus (Chironomidae), Culicoides sonorensis (Ceratapogonidae), and Polypedilum vanderplanki(Chironomidae), all of which have 3rd position GC content of 35% or less. We will replace them with protein sequences from Forpomyia taiwanana - GCA_963930915.1 (Ceratapogonidae), Belgica antarctica - GCA_000775305.1 (Chironomidae), and a transcriptome from a Similium sp. - GGBP00000000 (Simuliidae). All three have third codon positon GC sequences more similar to mosquitoes than to the other outgroups (A below).

Phylogenetic inference with IQTree on the new alignments including these outgroups results in essentially the same topology (D, below) as amino acids, NT1/2, or recoding NT3(C, below):

# Ultraconserved Elements
While analysis of protein coding genes from genomes or transcriptomes is quite common, one set of alternative markers that are regularly used are [ultraconserved elements](https://www.ultraconserved.org/). These highly-conserved areas are often analyzed with flanking sequence data. Previously analyses have shown that many UCE sets are located predominantly in protein-coding regions (>80% in flies - See [Van Dam et al. 2021](https://academic.oup.com/sysbio/article/70/2/307/5880562)). However, UCEs are difficult to partition - often, flanking regions include both coding and non-coding sequences, and codon positions are difficult to ascertain without additional steps, even if one knows what UCEs are protein coding and not. Moreover, non-coding regions may be subject to even more saturation than coding regions, in the same way that NT3 is typically more saturated than other codon positions.  None the less, UCEs allow for consistency recovery of phylogenetically informative loci across a wide range of species, and are substantially cheaper to generate than whole genomes or transcriptomes. Alternatives to UCEs with similar cost include hybrid or sequence capture and anchor hybrid enrichment, which often specifically target protein coding regions.

As the usage of UCEs is extremely common in many insects - and they were used by Pierce et al. 2025 - we wished to evaluate UCEs from mosquito genomes, see what phylogeny resulted from their analysis, and evaluate if they, too, showed issues with compositional heterogeineity. To do this, we used Dataseta C, described above, and augmented it with relevant outgroups for this analysis (Culicimorph outgroups).

## Extracting UCEs from Mosquito Genomes

## Determining if a UCE is Coding or Not

## BUSCO and UCE Overlap
As the vast majority of UCEs are in coding regions, one might wonder... are they overlapping with BUSCO genes?


## Alignment and Phylogenetic Analyses of UCEs
Once we had our UCEs extracted from genomes and the relevant categories, we started by aligning every UCE (mafft --genafpair --maxiterate 1000 --adjustdirection) and trimmed them with trimal... Following this, we further partitioned our UCEs that had coding regions into thirds. 

Although we don't know which nucleotide position corresponds to which codon position, this allows for us to consider each position in each UCE separately, such as for substitution model estimation, or to evaluate it based on something like GC content.

When we estimate the phylogeny from Dataset C, partitioned by nucleotide position within coding UCEs, we do recover a topology that has the Culicinae as non-monophyletic:

Again the question is... Why?

## Compositional Heterogeneity in UCEs

If you've read this far, you probably can already guess where this is going. For each taxa at each UCE (and at each position for coding UCEs), we calculated the GC content, and then we also summarized GC content for each UCE (and at each position for coding UCEs). As there is obvious differences in GC content between codon positions (see above), we then ranked nucleotide positions based on GC content, such that if in uce1, the first nucleotide position had a GC content of 35%, the second had a GC content of 60%, and the third had a GC content of 48%, we would rank the second nucleotide as "1", the third nucleotide as "2", and the first nucleotide as "3". We don't know which nucleotide position corresponds to which codon position, but given previous results that GC content in mosquito genes varies with codon position, we could rank all three 

We can do this with a handy R script - this assumes you have imported a dataframe that has a column that defines the partition called partition and defines the partition as marker_position, a Sequence_ID column, and a GC column. It will then output a "Rank" column for any partition with 1, 2, or 3, or "non" for others.
```
library(dplyr)
library(stringr)
df_processed <- df %>%
  mutate(
    name = str_extract(partition, "^[^_]+"),
    position = ifelse(str_detect(partition, "_"),
                      str_extract(partition, "(?<=_)\\d+"),
                      "non")
  ) %>%
  mutate(position = ifelse(position == "non", position, as.character(position))) %>%
  group_by(name, position) %>%
  mutate(avg_gc = mean(GC, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(name) %>%
  mutate(
    rank = ifelse(position == "non",
                  "non",
                  dense_rank(desc(avg_gc))))%>%
  ungroup() %>%
  select(-avg_gc)
```

Now this script automatically creates a column called "rank" which either includes the numeric rank from one to three where one has the highest GC content within that gene and three has the lowest GC content, or "non" for non-coding UCEs (as they are not split into nucleotide positions). We can then plot these data as follows:

```
library(ggplot2)
ggplot(df_processed,aes(x=GC,y=factor(Sequence_ID),fill=as.factor(rank)))+
  geom_boxplot(show.legend=FALSE,outliers=F,median.color = "white")+
  scale_fill_viridis(discrete=TRUE)+
  theme_bw()+ylab("Taxa")+
  xlab("GC")+facet_wrap(vars(factor(rank, level=rank_order)),nrow=1)
```

Please note that we can not directly interpret these as codon positions and so I do not suggest doing so. Indeed, this is a rough estimate of the idea there would be three "ranks" of composition within each coding gene. This is not always true, and for species where GC usage is not as distinct as mosquitoes, this will not necessarily be reliable. Though, it is kinda neat to see intron GC usage is so high... Makes one wonder... Anyway!

## Recoding UCEs

Given the compositional heterogeneity described above, we could recode UCEs as described in above sections.  What is the result when we do so, and then run a similar phylogenetic analysis on these recoded data as elsewhere?
<TWO SH>
The two subfamilies are once again represented in the data.  RY effectively mitigates compositional heterogeneity. 

# Concluding remarks
Compositional heterogeneity is generally recognized as a potential major source of error in phylogenetic inference. In the case of mosquitoes, the significant variability in GC usage at third codon position creates long-branch attraction artifacts when analyzing this position alone. While the results above demonstrate why NT3 can be problematic, it is worth noting it is none the less a useful source of information among closely related species.
