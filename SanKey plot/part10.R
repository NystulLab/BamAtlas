library(networkD3)
library(stringr)
library(dplyr)

MisLinks <- read.table("all.stat.xls", header = T, stringsAsFactors = F,sep="\t")
lev1 <- str_sort(unique(MisLinks$source),numeric = TRUE)
lev2 <- str_sort(unique(MisLinks$target),numeric = TRUE)
a <- setdiff(lev1,lev2)
b <- intersect(lev1,lev2)
c <- setdiff(lev2,lev1)

MisNodes <- data.frame(name = c(a,b,c))
MisNodes$group_color <- c(rainbow(length(a)),rainbow(length(b)),rainbow(length(c)))
Node2index = list()
Node2index[as.vector(MisNodes$name)] = 0:length(MisNodes$name)

MisLinks = MisLinks %>%
  mutate(source2 = unlist(Node2index)[as.vector(source)]) %>%
  mutate(target2 = unlist(Node2index)[as.vector(target)])

MisLinks$source <- factor(MisLinks$source,levels = c(a,b))
MisLinks$target <- factor(MisLinks$target,levels = c(b,c))
MisLinks <- arrange(MisLinks,source,target)
color <- MisNodes$group_color
names(color) <- MisNodes$name
MisLinks$group_color2 <- color[MisLinks$source] 

# painting
p <- sankeyNetwork(Links = MisLinks, 
              Nodes = MisNodes,
              Source = "source2", 
              Target = "target2",
              Value ="value",
              NodeID = "name",
              NodeGroup = "group_color", 
              LinkGroup  = "group_color2", 
              fontSize = 10
)

# saving
html_out <- paste0(out_dir, "/sankey.html")
pdf_out  <- paste0(out_dir, "/sankey.pdf")

htmlwidgets::saveWidget(p, file = html_out, selfcontained = F, libdir = paste0(out_dir,"/lib"))
webshot::webshot(url = html_out, file = pdf_out)
