### Define input ####

# input files
library(ape)
library(phytools)
library(stringr)
library(data.table)
library(dplyr)
library(glue)
f = glue

prefix    = "results/phylo/"
phy_list  = c("results/phylo/coeae1f_focal.fasta.treefile", "results/phylo/coeae1f_upstream.fasta.treefile", "results/phylo/coeae1f_downstream.fasta.treefile")
phy_name  = c("Focal haplotype","Upstream","Downstream")


clusters_df = fread("results/haplotype_clusters_metadata_r.tsv")
palette = clusters_df[,c('colour', 'Haplotype cluster')]

name='coeae1f'
region='focal'
phi = f("results/phylo/{name}_{region}.fasta.treefile")
phy = read.tree(phi)
phy_p = midpoint.root(phy)

meta_path = f("results/phylo/{name}_{region}.metadata.tsv")
meta = fread(meta_path) %>% as.data.frame()
meta = meta %>% arrange(factor(hap, levels = phy_p$tip.label))


#table(meta$karyotype, meta$`Sweep IDs`)

#meta_sweep3 = meta %>% filter(`Sweep IDs` == 3)
#meta_sweep3 %>% fwrite(., "results/sweep3_meta.tsv", sep="\t")


meta %>% with(which(aim_species == 'mela'))




plot_phylo = function(phi, meta, var, region, name, min_edge_length = 5e-5){
  
  phi = f("results/phylo/{name}_{region}.fasta.treefile")
  meta_path = f("results/phylo/{name}_{region}.metadata.tsv")
  
  # tree to distance matrix
  phy    = read.tree(phi)
  dis    = cophenetic.phylo(phy)
  meta = fread(meta_path) %>% as.data.frame()
  phy_p = midpoint.root(phy)
  meta = meta %>% arrange(factor(hap, levels = phy_p$tip.label))
  
  #phy = root(phy, outgroup=3368)
  
  phy_p$tip_label_sep  = gsub("_"," ",phy_p$tip.label)
  phy_p$tip_species    = meta$aim_species
  phy_p$tip_country = meta$country
  phy_p$tip_cluster   = ""
  phy_p$tip_karyotype  = ""
  #phy_p$tip_color     = brewe
  phy_p$tip.label      = rep("·", length(phy$tip.label))
  #phy_p$edge.length[phy_p$edge.length >  0.005] = 0.005
  phy_p$edge.length[phy_p$edge.length == 0]    = min_edge_length    #5e-5
  
  if (var == 'karyotype'){
    meta = meta %>% mutate("tip_color" = case_when(karyotype == '2l+a' ~ 'bisque2',
                                                   karyotype == '2la' ~ 'darkslategrey', TRUE ~ 'grey'))
  } else if (var == 'hap_cluster'){
    meta = meta %>% mutate("tip_color" = case_when(hap_cluster == 'C1' ~ '#7fffd4',
                                                 hap_cluster == 'C2' ~ '#000080',
                                                 hap_cluster == 'C3' ~ '#4169e1',
                                                 hap_cluster == 'C4' ~ '#87ceeb',
                                                 hap_cluster == 'C5' ~ '#ff8c00',
                                                 hap_cluster == 'C6' ~ '#7fffd4',
                                                  hap_cluster == 'WT' ~ '#d3d3d3',
                                                 hap_cluster == "" ~ 'black',
                                                 hap_cluster %in% c("mela", "meru", "quad") ~ 'black',
                                                 TRUE ~ "grey"))
  } else if (var == 'aim_species'){
    meta = meta %>% mutate("tip_color" = case_when(aim_species == 'gambiae' ~ 'dodgerblue',
                                                   aim_species == 'coluzzii' ~ 'indianred',
                                                   aim_species == 'arabiensis' ~ 'aquamarine3',
                                                   aim_species == 'mela' ~ 'cornsilk',
                                                   aim_species == 'meru' ~ 'cornsilk4',
                                                   aim_species == 'quad' ~ 'darkolivegreen',
                                                   TRUE ~ 'grey'))
  }
  plot.phylo(phy_p, type="unr",
             use.edge.length=T, show.tip.label=T, show.node.label=F,
             tip.color = meta$tip_color, lab4ut = "axial",
             edge.color = "slategray3",
             font = 1, edge.width = 2, node.depth=1, cex=5,
             main=f("{name} {region} | {var}"), no.margin = F)
}


jpeg(file=f("results/phylo/{name}_karyo.phylo.jpg"), width=10, height=20)



plot_phylo(phi, meta, var='hap_cluster', region="focal", name=name)


dev.off()



par(mar=c(1,1,1,1))


for (region in c("focal", "downstream", "upstream")){
  for (var in c("aim_species", "karyotype", "hap_cluster")){
    png(file=f("results/phylo/figures/{name}_{var}_{region}.phylo.png"), width=10, height=10)
    plot_phylo(phi, meta, var=var, region=region, name=name)
    dev.off()
  }

}

?svg
unique(meta$cluster_id)

plot_phylo(phi, meta, var='cluster_id', region="focal", name=name)




plot_phylo(phi, meta, var='karyotype', region="downstream", name=name)
plot_phylo(phi, meta, var='aim_species', region="downstream", name=name)
plot_phylo(phi, meta, var='Sweep IDs', region="downstream", name=name)

plot_phylo(phi, meta, var='karyotype', region="upstream", name=name)
plot_phylo(phi, meta, var='aim_species', region="upstream", name=name)
plot_phylo(phi, meta, var='Sweep IDs', region="upstream", name=name)

# 6 sweeps
# 1 arabiensis - CNV
# 3 gambiae one of which has CNV
# 1 coluzzii
# 1 mixed gambiae and coluzzii


