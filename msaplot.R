library(tidyverse)

df=read_csv("mammals_table.csv")

library(rentrez)

#TODO: Add Entrez API key
ENTREZ_API_KEY=""

if (ENTREZ_API_KEY != "") {
  set_entrez_key(ENTREZ_API_KEY)
}

# Get protein sequences
ef=entrez_fetch("protein", id=df$protein_id, rettype = "fasta")
write(ef, "ace2.fasta")

library(Biostrings)
aa=readAAStringSet("ace2.fasta")

# Abbreviated naming
names(aa)=str_c(df$abbreviated, " (", df$common_name, ")")

library(msa)

# MSA
m=msa(aa, order = "input")
mm=msaConvert(m, type="bios2mds::align")
mm=lapply(mm, function(x) c(x[24:43], "-", x[81:83], "-", x[353:359]))
class(mm)="align"
bios2mds::export.fasta(mm, "ace2.aligned.fasta")

library(ggtree)

# Create tree
#d=as.dist(stringDist(subseq(aa, start = 1, end = min(width(aa))), method = "hamming")/width(aa)[1])
#t=ggtree(ape::bionj(d)) + geom_tiplab()

# Imported tree from timetree.org
dt=read.tree("mammals.nwk")

# Same naming scheme as MSA
dt$tip.label=replace(dt$tip.label,dt$tip.label == "Canis_lupus","Canis lupus familiaris")
dt$tip.label=str_replace(dt$tip.label, "_", " ")
dt$tip.label=unlist(lapply(dt$tip.label, function(x) str_c(df[df$scientific_name == x,]$abbreviated, " (", df[df$scientific_name == x,]$common_name, ")")))

# Scale down branches
#dt=treeio::rescale_tree(dt, 0.001)
dt$edge.length=dt$edge.length/100

t=ggtree(dt) 
#+ geom_tiplab(align = T, hjust = 1,offset = 25, linesize = 0.01)

# Put humans first
t=flip(t,27,17)
t=flip(t, 18,21)
t=flip(t, 1, 2)
#t=flip(t, 19, 20)
#tt=t + geom_tiplab(align=T) + xlim(0,2)


library(awtools)

mp=msaplot(t, "ace2.aligned.fasta", window = c(1,30), width = 2) + geom_tiplab(offset = 2.05) + xlim(0,5) + labs(fill="Amino acid") + geom_label2(aes(subset=(node==29), label="Primates"), nudge_x = -0.4) + geom_label2(aes(subset=(node==28), label="Rodents"), nudge_x = -0.3) + geom_label2(aes(subset=(node==21), label="Bats"), nudge_x = -0.1) + geom_label2(aes(subset=(node==20), label="Ungulates"))

mp=mp + geom_cladelabel(27, label = "chrX", offset = 3.6) + geom_cladelabel(20, label = "chrX", offset = 3.6) + geom_cladelabel(10, label = "chr9", offset = 3.6) + geom_cladelabel(6, label = "chr1", offset = 3.6) + geom_cladelabel(1, label = "chrX", offset = 3.6)

mp=mp + theme_tree2() + scale_x_continuous(breaks = c(1.7, 2.47, 2.8, 4.63), labels = c("24-43", "81-83", "353-357", "chromosome")) +
  scale_fill_manual(values=c("#ffffff",bpalette,"#ff7300"), labels=c("-","A","D","E","F","K","L","N","Q","R","S","Y","G","H","I","T","V","M"))

ggsave("ace2_homology.png", plot=mp, device = "png", dpi = 320, width=12, height = 10)