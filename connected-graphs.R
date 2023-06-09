corr_graph <- function(
  df,
  layout = "layout_nicely",
  output = c("windows", "quartz", "bmp", "jpeg", "png"),
  filename = "corr_graph",
  domains = NULL,
  res = 150,
  height = 6,
  width = 12,
  r.cut.off = 0.5,
  alpha.cut.off = 0.05,
  use.significance = FALSE,
  set.seed = 19062020,
  max.line.width = 6,
  vertex.cex = 0,
  label.cex = 0.65,
  border = "light grey",
  line.col = c("#3366FF", "#FA8072"),
  margin = 0.1,
  ncol = 2,
  nrow = 1,
  asp = 1,
  omit.rest = TRUE,
  use.only = NULL,
  scale.vertex = FALSE,
  legend.pos = c("none", "bottomleft", "bottomright", "topleft", "topright"),
  annotate = FALSE,
  send.to.output = TRUE) {
  
  # Note: see https://rstudio-pubs-static.s3.amazonaws.com/362044_903076131972463e8fdfcc00885fc9a6.html
  # to have a look at the available layouts in igraph (towards the bottom)

  try.pkg <- try(require(stats))
  if(!try.pkg) {
    install.packages("stats")
    require(stats)
  }
  try.pkg <- try(require(Hmisc))
  if(!try.pkg) {
    install.packages("Hmisc")
    require(Hmisc)
  }  
  try.pkg <- try(require(igraph))
  if(!try.pkg) {
    install.packages("igraph")
    require(igraph)
  }
  try.pkg <- try(require(scales))
  if(!try.pkg) {
    install.packages("scales")
    require(scales)
  }
  
  if(!is.null(use.only)) domains <- names(df) %in% use.only
  set.seed(set.seed)
  
  # Reordering the data produces a better circular and grid layout
  
  k <- hclust(as.dist(sqrt(2*(1-cor(na.omit(df))))), method = "centroid")$order
  df <- df[,k]
  if(!is.null(domains)) domains <- domains[k]
  
  # Calculate the correlations and use their absolute values
  # to create the graph / network
  
  correlation_matrix <- rcorr(as.matrix(df))
  
  M <- abs(correlation_matrix$r)
  diag(M) <- 0
  total_rho <- rowSums(M)
  if(!scale.vertex) total_rho <- rep(50, length(total_rho)) 
  
  if(!is.null(use.only)) {
    kk <- which(names(df) %in% use.only)
    M[-kk,-kk] <- 0
  }
  
  p <- correlation_matrix$P
  diag(p) <- 0

  ifelse(use.significance, M[p > alpha.cut.off] <- 0, M[M < r.cut.off] <- 0)
  
  k <- which(apply(M, 2, function(x) {
    all(x == 0)
  }))
  if(omit.rest & length(k)) {
    M <- M[-k,-k]
    total_rho <- total_rho[-k]
    if(!is.null(domains)) domains <- domains[-k]
  }
  
  M <- sqrt(M)
  net <- graph_from_adjacency_matrix(M, weighted=T, mode="undirected", diag=F)
  
  # Below is to retrieve whether the correlations are positive
  # or negative
  
  M <- correlation_matrix$r
  diag(M) <- 0
  
  if(!is.null(use.only)) M[-kk,-kk] <- 0
  
  ifelse(use.significance, M[p > alpha.cut.off] <- 0, M[abs(M) < r.cut.off] <- 0)

  if(omit.rest & length(k)) M <- M[-k,-k]
  M <- M[lower.tri(M)]
  M <- M[M != 0]
  
  # Set up the output device
  
  output <- output[1]
  ifelse(output == "windows" | output == "quartz",
    device <- paste0(output,"(height = ",height,", width = ",width,")"),
    device <- paste0(output,"(\"",filename,".",output,"\", height = ",
                     height,", width = ",width,", res = ",res,
                     ", units = \"in\")"))
  if(send.to.output) eval(parse(text = device))
  
  # Set the attributes
  
  l <- paste0(layout,"(net)")
  l <- eval(parse(text = l))

  if(!is.null(domains)) {
    ifelse(is.null(use.only), gp <- factor(domains),
                              gp <- factor(domains, levels = c("TRUE", "FALSE")))
    gp <- as.numeric(gp)
    try.pkg <- try(require(RColorBrewer))
    if(!try.pkg) {
      install.packages("RColorBrewer")
      require(RColorBrewer)
    }
    ifelse(is.null(use.only), bd.pal <- brewer.pal(max(gp), "Dark2"),
                              bd.pal <- c("black", NA))
    border <- bd.pal[gp]
  }
  
  
  V(net)$frame.color <- border
  V(net)$color <- scales::alpha("grey", 0.4)
  V(net)$size <- vertex.cex*total_rho
  V(net)$label.cex <- label.cex
  V(net)$label.color <- "black"
  
  pal <- ifelse(M > 0, line.col[1], NA)
  E(net)$color <- scales::alpha(pal, E(net)$weight)
  E(net)$curved <- 0.3
  w <- E(net)$weight^2
  w <- (max.line.width - 1) * (w - min(w)) / (max(w) - min(w)) + 1
  E(net)$width <- w
  corr <- paste0("+",round(E(net)$weight, 3))
  corr[M < 0] <- NA
  if(annotate) E(net)$label <- corr
  E(net)$label.cex <- 1.2*label.cex
  E(net)$label.color <- "black"

  
  # Plot the graphs
  
  par(mfrow=c(nrow,ncol))
  par(mai=rep(margin,4))
  
  plot(net, 
       layout = l,
       rescale= FALSE,
       xlim = c(min(l[,1]), max(l[,1])),
       ylim = c(min(l[,2]), max(l[,2])),
       asp = asp
  )
  
  if(!is.null(domains) & grepl("left", legend.pos[1])) {
    legend(legend.pos, legend = levels(as.factor(domains)),
           col = bd.pal, cex = 0.8, bty = "n", pch = 1)
  }
  
  pal <- ifelse(M < 0, line.col[2], NA)
  E(net)$color <- scales::alpha(pal, E(net)$weight)
  corr <- paste0("-",round(E(net)$weight, 3))
  corr[M > 0] <- NA
  if(annotate) E(net)$label <- corr
  
  plot(net, 
       layout = l,
       rescale= FALSE,
       xlim = c(min(l[,1]), max(l[,1])),
       ylim = c(min(l[,2]), max(l[,2])),
       asp = asp
  )
  
  if(!is.null(domains) & grepl("right", legend.pos[1])) {
    legend(legend.pos, legend = levels(as.factor(domains)),
           col = bd.pal, cex = 0.8, bty = "n", pch = 1)
  }
  
  if(output != "windows" & output != "quartz" & send.to.output) dev.off()
  
}