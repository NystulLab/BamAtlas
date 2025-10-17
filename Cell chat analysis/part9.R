library(CellChat)

dt <- netVisual_bubble(cellchat, comparison = c(1, 2), return.data = TRUE)
df <- dt[['communication']][, c(-5,-6,-12)]

Signaling <- head((df %>% group_by(pathway_name) %>% summarise(prob = sum(prob.original)) %>% arrange(-prob))$pathway_name, 5)

color.text <- c("#377EB8","#FF7F00")
p <- netVisual_bubble(cellchat, comparison = c(1, 2), angle.x = 45, signaling = Signaling, remove.isolate = T, color.text = color.text, return.data = F)

ggsave(p, file = "Diff.CommunProb.top5.bubble.pdf", width = 7, height = 7)
