

rm(list=ls())

library(dplyr)
library(ggplot2)
library(plotly)
library(scales)
library(RColorBrewer)
library(patchwork)
library(reshape)
library(sm)
library(tidyr)

dir.create( path = paste0( getwd(), '/plots/data/diet'), showWarnings = TRUE, recursive = TRUE)

# Diet data ----------------------

diet <- read.csv2( './data/Diet/hke_completo.csv')[,-c(1,5)]

colnames(diet) <- c('Prey','Size','Percentage','Year')

diet$Percentage <- as.numeric(diet$Percentage)
diet$Size <- factor(diet$Size, levels = c('< 18', '18 - 34', '> 34'))
diet$Prey <- as.factor(diet$Prey)

diet_complete <- diet %>% complete( Prey, Size, Year, fill = list(Percentage = 0))

yplot <- ggplot( subset(diet_complete, Year>2008), aes(x = Size, y = Percentage, fill = Prey)) +
  facet_wrap(~ Year, ncol = 5) +
  geom_bar(stat = "identity", position = "stack") + labs(title = "Percentage by Prey and Size") +
  theme_minimal()
yplot


# By Prey and Size -------------------

by_prey_size <- diet_complete %>%
  group_by( Prey, Size) %>%
  summarise( Percentage = mean(Percentage, na.rm = TRUE), .groups = 'drop')

head(by_prey_size)

gplot <- ggplot( by_prey_size, aes(x = Size, y = Percentage, fill = Prey)) +
  geom_bar(stat = "identity", position = "stack") + 
  labs(title = "Percentage by Prey and Size") +
  theme_minimal()
gplot



# By Prey only -----------------------------------

by_prey <- by_prey_size %>%
  group_by( Prey) %>%
  summarise( Mean_Percentage = mean(Percentage, na.rm = TRUE), .groups = 'drop')

qplot <- ggplot( by_prey, aes(x = "", y = Mean_Percentage, fill = Prey)) +
  geom_bar(stat = "identity", width = 1) + coord_polar(theta = "y") +
  labs(title = "Percentage by Prey") + theme_minimal() + 
  theme( axis.title.x = element_blank(), axis.title.y = element_blank(), 
         axis.text = element_blank(),panel.grid = element_blank())
qplot



# Adding 'others' ----------------------------

top_prey <- diet_complete %>%
  group_by(Prey) %>%
  summarise(Total_Percentage = sum(Percentage, na.rm = TRUE)) %>%
  arrange(desc(Total_Percentage)) %>%
  slice_head(n = 11) %>%
  pull(Prey)  

diet_complete2 <- diet_complete %>%
  mutate(Prey = ifelse(!Prey %in% top_prey | Prey == "otros", "Other", as.character(Prey))) %>%
  group_by(Prey, Size, Year) %>%
  summarise(Percentage = sum(Percentage, na.rm = TRUE), .groups = 'drop')

yplot2 <- ggplot( subset(diet_complete2, Year>2008), aes(x = Size, y = Percentage, fill = Prey)) +
  facet_wrap(~ Year, ncol = 5) +
  geom_bar(stat = "identity", position = "stack") + labs(title = "Percentage by Prey and Size") +
  theme_minimal()
yplot2

by_prey_size2 <- diet_complete %>%
  group_by(Prey, Size) %>%
  summarise(Percentage = mean(Percentage, na.rm = TRUE), .groups = 'drop') %>%
  mutate(Prey = ifelse(Percentage < 0.025 | Prey == "otros", "Other", as.character(Prey))) %>%
  group_by(Prey, Size) %>%
  summarise(Percentage = sum(Percentage, na.rm = TRUE), .groups = 'drop')

gplot2 <- ggplot( by_prey_size2, aes(x = Size, y = Percentage, fill = Prey)) +
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = brewer.pal(n = 12, "Set3")) + 
  labs(title = "Percentage by Prey and Size") +
  theme_minimal()
gplot2

by_prey2 <- by_prey %>%
  mutate(Prey = ifelse(Mean_Percentage < 0.025 | Prey == "otros", "Other", as.character(Prey))) %>%
  group_by(Prey) %>%
  summarise(Mean_Percentage = sum(Mean_Percentage, na.rm = TRUE)) %>%
  arrange(desc(Mean_Percentage))
by_prey2

ggsave( "./plots/data/diet/diet1.jpg", width = 6, height = 4)

qplot2 <- ggplot(by_prey2, aes(x = "", y = Mean_Percentage, fill = Prey)) +
  geom_bar(stat = "identity", width = 1, color = "white") + 
  coord_polar(theta = "y") +
  scale_fill_manual(values = brewer.pal(n = min(11, length(unique(by_prey2$Prey))), "Set3")) + 
  labs(title = "Percentage by Prey") +
  theme_minimal() + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
                          axis.text = element_blank(),panel.grid = element_blank())
qplot2

ggsave( "./plots/data/diet/diet2.jpg", width = 6, height = 4)



# Cannibalism -------------------------------------------

cannibal_complete <- subset( diet_complete, Prey == 'M.merluccius')
cannibal_complete

cannibal_byyear <- cannibal_complete %>%
  group_by( Prey, Year) %>%
  summarise( Percentage = mean(Percentage, na.rm = TRUE), .groups = 'drop')

cannibal <- subset( by_prey_size, Prey == 'M.merluccius')
cannibal

cannibal2 <- subset(cannibal_complete, Year>=1996)

# cannibalplot <- ggplot( subset(cannibal_complete, Year>=1996), aes(x = Year, y = Percentage, fill = Size)) +
#   geom_bar(stat = "identity", position = "stack") + theme_minimal() +
#   geom_line( data = cannibal_byyear, x = Year, y = Percentage, linetype = 2) +
#   geom_hline( data=cannibal, aes( yintercept = Percentage, color = Size), linetype = 2) +
#   labs( title = "Cannibalism", fill = 'Predator Size', color = 'Average')
# 
# cannibalplot

cannibalplot <- ggplot( cannibal_complete, aes(x = Year, y = Percentage, fill = Size)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) + 
  theme_minimal() +
  geom_line(data = cannibal_byyear, aes(x = Year, y = Percentage), color = "black", linetype = 1, inherit.aes = FALSE) +
  geom_hline(data = cannibal, aes(yintercept = Percentage, color = Size), linetype = 2) +
  labs(title = "Cannibalism", fill = 'Predator Size', color = 'Average')

print(cannibalplot)

cannibal_comp2 <- cannibal_complete
cannibal_comp2$Percentage <- cannibal_complete$Percentage/3

cannibalplot <- ggplot( cannibal_comp2, aes(x = Year, y = Percentage, fill = Size)) +
  geom_bar(stat = "identity", position = 'stack') + 
  theme_minimal() +
  # geom_hline(data = cannibal, aes(yintercept = Percentage, color = Size), linetype = 2) +
  labs(title = "Cannibalism", fill = 'Predator Size', color = 'Average')

cannibalplot

ggsave( "./plots/data/diet/cannibalism.jpg", width = 6, height = 4)


# Save --------------

pdf("./plots/data/diet/diet.pdf", width = 10, height = 6, onefile = TRUE)
print(yplot2)
print(gplot2)
print(qplot2)
print(cannibalplot)
dev.off()

save( by_prey, by_prey_size, diet, diet_complete, by_prey2, by_prey_size2, 
      diet, diet_complete2, cannibal, cannibal_complete, cannibal_byyear,
      file = './data/Diet.RData')

