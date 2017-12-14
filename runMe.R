library("poppr")
library("tidyverse")
library("igraph")
library("ggraph")
library("visNetwork")
dat11 <- read_csv("https://osf.io/cmnk2/download")
loc   <- read_csv("https://osf.io/6huw3/download")$Locus
the_loci <- which(gsub("[(FH)]", "", names(dat11)) %in% loc)
dat11 <- df2genind(dat11[the_loci], ind.names = dat11$Isolate, strata = dat11[1:7], ploidy = 1)
dat11 <- as.genclone(dat11)
mll.custom(dat11) <- strata(dat11)$MCG
mcgmlg <- as.data.frame(table(mll(dat11, "original"), mll(dat11, "custom"))) %>%
  setNames(c("MLG", "MCG", "Freq")) %>%
  mutate(MLG = as.character(MLG)) %>%
  mutate(MCG = as.character(MCG)) %>%
  as_tibble() %>%
  filter(Freq > 0)

make_mcgmlg_graph <- function(x){
  gdf <- mutate(x, MLG = paste0('MLG.', MLG))
  MLGS <- gdf %>% 
    group_by(MLG) %>%
    summarize(size = sum(Freq)) %>%
    rename(vertex = MLG)
  MCGS <- gdf %>% 
    group_by(MCG) %>%
    summarize(size = sum(Freq)) %>%
    rename(vertex = MCG)
  VAT <- bind_rows(MLGS, MCGS)
  g <- gdf %>% 
    select(MLG, MCG, Freq) %>%
    rename(weight = Freq) %>%
    graph_from_data_frame(vertices = VAT)
  V(g)$type <- ifelse(grepl("MLG", V(g)$name), "Multilocus Genotype", "Mycelial Compatibility Group")
  g
}

add_communities <- function(g, clusters){
  comm <- communities(clusters)
  commlist <- setNames(rep(names(comm), lengths(comm)), unlist(comm, use.names = FALSE))
  V(g)$community <- commlist[V(g)$name]
  g
}
g <- make_mcgmlg_graph(mcgmlg)
g <- add_communities(g, cluster_walktrap(g))
osize <- V(g)$size
set.seed(2017-05-03)
# set.seed(2017-08-02)

lay2 <- create_layout(g, layout = "igraph", 
                      algorithm = "fr", 
                      weights = rep(1, length(E(g)$weight))
                      # weights = ifelse(E(g)$weight == 1, 1, 1 + E(g)$weight/100)
)
the_communities <- data_frame(vertex = lay2$name, community = lay2$community) %>%
  mutate(comm = case_when(
    .$community == 7 ~ "A",
    .$community == 5 ~ "B",
    .$community == 1 ~ "C",
    TRUE ~ as.character(.$community)
  )) %>%
  group_by(community) %>%
  mutate(comm = ifelse(n() > 10, paste("Community", comm), "Other Communities (n < 10)")) %>%
  mutate(comm = comm) 


make_rgb <- function(x, alpha = 1){
  out <- col2rgb(x, alpha = TRUE)
  out[4, ] <- floor(out[4, ]*alpha)
  out      <- apply(out, 2, paste, collapse = ", ")
  paste0("rgba(", out, ")")
}
community_colors <- viridis::viridis(4, direction = -1)[as.integer(factor(the_communities$comm))]
vg <- g %>%
  set_vertex_attr("size", value = osize * 10) %>%
  set_vertex_attr("value", value = osize) %>%
  set_vertex_attr("label", value = NULL) %>%
  set_vertex_attr("color", value = community_colors) %>%
  set_vertex_attr("shape", value = ifelse(V(.)$type == "Multilocus Genotype", "triangle", "dot")) %>%
  set_edge_attr("width", value = E(.)$weight) %>%
  toVisNetworkData()
vg$nodes <- vg$nodes %>% 
  dplyr::group_by(id) %>%
  dplyr::mutate(color = list(list(background = make_rgb(color, 0.8), 
                                  border = make_rgb(rep("black", n()), 0.8),
                                  highlight = list(background = make_rgb(color),
                                                   border = make_rgb(rep("black", n()))
                                  )
  ))) %>%
  dplyr::mutate(title = paste0("<p>",
                               "<i>", type, " ", gsub("MLG.", "", id), "</i><br>",
                               "<b>Community: ", community, "</b><br>",
                               "<b>N: ", size/10, "</b><br>",
                               "</p>")
  )
vg$edges <- vg$edges %>%
  dplyr::mutate(title = paste0("<p>",
                               "<b>Mycelial Compatibility Group: ", from, "</b><br>",
                               "<b>Multilocus Genotype: ", gsub("MLG.", "", to), "</b><br>",
                               "<b>N: ", weight, "</b>",
                               "</p>")
  )
vgn <- visNetwork(nodes = vg$nodes, edges = vg$edges,# height = "500px", 
                  main = "Relation of Multilocus Genotypes and MCGs")
set.seed(2017-05-03)

vgn %>%
  # visIgraphLayout("layout_nicely") %>% # activate this for a poseable network
  visNodes(label = NULL, shadow = TRUE) %>%
  visOptions(highlightNearest = list(enabled = TRUE, hover = TRUE), 
             nodesIdSelection = TRUE)
