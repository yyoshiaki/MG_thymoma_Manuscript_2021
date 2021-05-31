library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggdendro)
library(cowplot)
library(patchwork) 
library(reshape2)
library(RColorBrewer)
library(ggnewscale)
citation("ggnewscale")
setwd("~/bioinformatics/Train/MG_thymoma/")


data <- read.table("./minor_cluster_module_stats.csv", sep = ",", 
                   header = TRUE, row.names = NULL, stringsAsFactors =TRUE)
colnames(data)
names(data)[1] <- "clusters"

color <- c("black", "turquoise", "blue", "yellow", "green", "red", "grey", 
           "autoantibodies", "GWAS")
means <- paste("mean_", color, sep = "")
padjs <- paste("padj_", color, sep = "")

mean.data <- data[, c("clusters", means, padjs)]
mean.data

clusters <- factor(c(row.names(mean.data)))

# convert to tidydata
long_data <- mean.data %>% gather(means, key = "modules", value = "exp") %>% 
  gather(padjs, key = "padjs", value = "padj")
# before using abs(), adding categorical data which negative or positive sign
long_data <- long_data %>% 
  mutate(sign = case_when(exp >= 0 ~ 1, exp < 0 ~ -1))
  
# convert zero padjs to 1e-50  
long_data$padj <- ifelse(long_data$padj== 0, 1e-8, long_data$padj)

# adding -log10 columns
long_data$logs <- -(log10(abs(long_data$padj))) 

long_data$log_padj <- (long_data$logs * long_data$sign)

long_data$modules <- factor(long_data$modules, levels=unique(long_data$modules))
long_data$clusters <- factor(long_data$clusters, levels=(unique(long_data$clusters)))
levels(long_data$clusters)

# pdf("./210510_dotplot_3.pdf", width = 4, height = 8)

# long_data %>% 
#   group_by(clusters) %>% 
#   summarise(amt = sum(exp)) 

# value <- seq(long_data$log_padjs)

pdf("./210531_dotplot_7.pdf", width = 4.5, height = 8)
p <- ggplot(data = long_data, aes(x = modules, y = clusters, colour = log_padj)) + 
  geom_point(aes(size = abs(exp)), shape = 20, alpha = 0.5) + 
  # scale_color_manual(values =  c("blue", "red")) +
  scale_colour_gradientn(colours = c("blue","white","red"), limits = c(-8, 8)) +
  scale_size(range = range(c(-10:16))) + 
  theme_bw()

p

dev.off()

# scale_colour_gradient2(low = "blue", high = "red", na.value = NA, aesthetics = "colour") + 


# p <- ggplot() + 
#   geom_point(data = negs, aes(x = modules, y = clusters, size = logs, color = logs, fill = logs), shape = 21, alpha = 1) + 
#   scale_colour_gradient2(low = "white", high = "blue", na.value = NA, aesthetics = "fill") + 
#   scale_size(range = c(-8, 8))  +
#   new_scale_color() +
#   new_scale_fill() + 
#   geom_point(data = pos, aes(x = modules, y = clusters, size = logs, color = logs, fill = logs), shape = 21, alpha = 1) + 
#   scale_color_gradient2(low = "white", high = "red", na.value = NA, aesthetics = "fill") 

p

# dev.off()

# long_data$modules <- as.character(long_data$modules)

# mean.data <- melt(mean.data, id = modules)
# mean.data <- transform(mean.data, clusters = data$clusters)

# modules_ <- factor(x = modules, levels = color)
# clusters_ <- factor(x = data$clusters, levels = data$clusters)

# long_data$clusters <- factor(long_data$clusters, 
# levels = rev(levels(long_data$clusters)))

# neg <- (mean_long$exp < 0)
# negs <- mean_long[neg,]
# negs[,"sign"] = as.factor(c("neg"))
# negs$logs <- -(log10(abs(negs$mean_exp)))
# nrow(negs)
# neg <- as.data.frame(neg)
# 
# pos <- (mean_long$mean_exp >= 0)
# pos <- mean_long[pos,]
# pos[,"sign"] = as.factor(c("pos"))
# pos$logs <- -(log10(abs(pos$mean_exp)))
# nrow(pos)
# pos <- as.data.frame(pos)
# 
# means <- rbind(negs, pos)

# p <- ggplot() + 
#   geom_point(data = negs, aes(x = modules, y = clusters, size = logs, color = logs, fill = logs), shape = 21, alpha = 1) + 
#   scale_colour_gradient2(low = "white", high = "blue", na.value = NA, aesthetics = "fill") + 
#   scale_size(range = c(-15, 15))  +
#   new_scale_color() +
#   new_scale_fill() + 
#   geom_point(data = pos, aes(x = modules, y = clusters, size = logs, color = logs, fill = logs), shape = 21, alpha = 1) + 
#   scale_color_gradient2(low = "white", high = "red", na.value = NA, aesthetics = "fill") 
# p
