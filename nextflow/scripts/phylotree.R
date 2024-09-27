library(treedataverse)
library(ape)
library(ggtree)
mytree_newick <- ape::read.tree("/diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/intermediate/SB28_DAISY_DAY0_sub2_A3_S122_L007/BC2.newick")

pdf("/diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/results/BC2_tree_newick.pdf", width = 12, height = 7)
plot(mytree_newick)
dev.off()


mytree_iqtree <- read.tree("/diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/intermediate/SB28_DAISY_DAY0_sub2_A3_S122_L007/BCall_test.nexus.treefile")

p <- ggtree(mytree_iqtree, layout = 'circular', branch.length = 'none')

pdf("/diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/results/all_tree_test.pdf", width = 50, height = 50)
p
dev.off()

### Testing for circular
data(iris)
rn <- paste0(iris[,5], "_", 1:150)
rownames(iris) <- rn
d_iris <- dist(iris[,-5], method="man")

tree_iris <- ape::bionj(d_iris)
grp <- list(setosa     = rn[1:50],
            versicolor = rn[51:100],
            virginica  = rn[101:150])

p_iris <- ggtree(tree_iris, layout = 'circular', branch.length='none') + geom_tiplab() 
groupOTU(p_iris, grp, 'Species') + aes(color=Species) +
  theme(legend.position="right")

pdf("/diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/results/iris_testing.pdf", width = 12, height = 7)
p_iris
dev.off()

# Libraries
library(ggraph)
library(igraph)
library(tidyverse)
library(viridis)
 
edges <- flare$edges
vertices <- flare$vertices
mygraph <- graph_from_data_frame(edges, vertices=vertices)

p <- ggraph(mygraph, layout = 'circlepack', weight=size) + 
  geom_node_circle(aes(fill = depth))

pdf("/diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/results/edges_vertices_testing.pdf", width = 20, height =20)
p
dev.off()

# We need a data frame giving a hierarchical structure. Let's consider the flare dataset:
edges2 <- read.table("/diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/intermediate/edges_5.txt", header = FALSE, sep = "\t")
colnames(edges2) <- c("from","to")
# edges2$from <- paste("BC.",edges2$from,sep="")
# edges2$to <- paste("BC.",edges2$to,sep="")

vertices2 <- read.table("/diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/intermediate/vertices_5.txt", header = FALSE, sep = "\t")
colnames(vertices2) <- c("name","size")
vertices2 <- vertices2[-8117,]

g <- graph_from_data_frame(edges2, vertices=vertices2, directed = TRUE)
 
p2 <- ggraph(g, layout = 'circlepack', weight=size) + 
  geom_node_circle(aes(fill = depth))

# 500 header
# We need a data frame giving a hierarchical structure. Let's consider the flare dataset:
edges2 <- read.table("/diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/intermediate/edges_500_3.txt", header = FALSE, sep = "\t")
colnames(edges2) <- c("from","to")

vertices2 <- read.table("/diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/intermediate/static_BC_500_2.txt", header = FALSE, sep = "\t")
colnames(vertices2) <- c("name","size")

# vertices2$size <- vertices2$size +1 

g <- graph_from_data_frame(edges2, vertices=vertices2, directed = TRUE)
 
p2 <- ggraph(g, layout = 'circlepack', weight=size) + 
  geom_node_circle(aes(fill = size), alpha = 0.3) +
  scale_fill_gradientn(colours = c("white","red","red","red","red","red","red","red")) +
  # scale_fill_manual(values=c("0" = "white", "1" = viridis(4)[1], "2" = viridis(4)[2], "3" = viridis(4)[3], "4"=viridis(4)[4])) +
  # scale_color_manual( values=c("0" = "white", "1" = "black", "2" = "black", "3" = "black", "4"="black") ) +
  theme_void()

pdf("/diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/results/edges_vertices_alpha_2.pdf", width = 30, height =30)
p2
dev.off()