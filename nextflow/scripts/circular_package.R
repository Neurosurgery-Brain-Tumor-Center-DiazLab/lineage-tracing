### Making a circule packege
# awk -F "," '{print "BC."$3"\tBC."$3"."$4}' SB28_DAISY_DAY13_sub2_A4_S123_L007_merged_homdomain_10x_barcode_combined_max_10x_UMI_static_mut_clean_remove_PCR_umi_col.csv  | grep -vi "static" > SB28_DAISY_DAY13_sub2_A4_S123_L007_edges.txt
# cut -f1 SB28_DAISY_DAY13_sub2_A4_S123_L007_edges.txt | sort | uniq | awk '{print $1"\t0"}' > SB28_DAISY_DAY13_sub2_A4_S123_L007_vertices_append.txt
# awk '{print $2}' SB28_DAISY_DAY13_sub2_A4_S123_L007_edges.txt | sort | uniq -c | awk '{print $NF"\t"$1}' > SB28_DAISY_DAY13_sub2_A4_S123_L007_vertices.txt
# sort SB28_DAISY_DAY13_sub2_A4_S123_L007_edges.txt | uniq > SB28_DAISY_DAY13_sub2_A4_S123_L007_edges_2.txt
# cat SB28_DAISY_DAY13_sub2_A4_S123_L007_vertices.txt SB28_DAISY_DAY13_sub2_A4_S123_L007_vertices_append.txt > SB28_DAISY_DAY13_sub2_A4_S123_L007_vertices_2.txt
library(ggraph)
library(igraph)
library(tidyverse)
library(viridis)

edges2 <- read.table("/diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/intermediate_2/SB28_DAISY_DAY13_sub2_A4_S123_L007_edges_2.txt", header = FALSE, sep = "\t")
colnames(edges2) <- c("from", "to")
# edges2$from <- paste("BC.",edges2$from,sep="")
# edges2$to <- paste("BC.",edges2$to,sep="")

vertices2 <- read.table("/diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/intermediate_2/SB28_DAISY_DAY13_sub2_A4_S123_L007_vertices_2.txt", header = FALSE, sep = "\t")
colnames(vertices2) <- c("name", "size")

g <- graph_from_data_frame(edges2, vertices = vertices2, directed = TRUE)

p2 <- ggraph(g, layout = "circlepack", weight = size) +
  geom_node_circle(aes(fill = size), alpha = 0.5) +
  scale_fill_gradientn(colours = c(
    "white", "red", "red", "red", "red", "red", "red", "red", "red",
    "red", "red", "red", "red", "red", "red", "red", "red",
    "red", "red", "red", "red", "red", "red", "red", "red",
    "red", "red", "red", "red", "red", "red", "red", "red",
    "red", "red", "red", "red", "red", "red", "red", "red"
  )) +
  # scale_fill_manual(values=c("0" = "white", "1" = viridis(4)[1], "2" = viridis(4)[2], "3" = viridis(4)[3], "4"=viridis(4)[4])) +
  # scale_color_manual( values=c("0" = "white", "1" = "black", "2" = "black", "3" = "black", "4"="black") ) +
  theme_void()

pdf("/diazlab/data3/.abhinav/projects/DAISY/lineage-tracing/results/SB28_DAISY_DAY13_edges_vertices_alpha_2.pdf", width = 30, height = 30)
p2
dev.off()
