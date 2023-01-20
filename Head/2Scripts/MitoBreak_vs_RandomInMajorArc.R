library(ggpubr)


data <- read.csv("MitoBreakDB_121219_good.csv", sep=';')

xplot <- ggplot(data, aes(Middle)) +
  geom_histogram(alpha=0.4, position = 'dodge') + theme_minimal() +
  xlab('') + ylab('') + scale_x_continuous(breaks=c(6000, 9000, 12000, 15000)) +
  scale_y_continuous(breaks = c(100, 200, 300, 400, 500, 600)) +
  scale_fill_manual(values=c('red', 'grey'))

xplot



data2 <- read.csv("MitoBreakDB_121219_rand.csv", sep=';')


yplot <- ggplot(data2, aes(Middle)) +
  geom_histogram(alpha=0.4, position = 'dodge') + theme_minimal() +
  xlab('') + ylab('') + scale_x_continuous(breaks=c(6000, 9000, 12000, 15000)) +
  scale_y_continuous(breaks = c(100, 200, 300, 400, 500, 600)) +
  scale_fill_manual(values=c('red'))

# Cleaning the plots

# Arranging the plot
ggarrange(xplot, NULL, yplot, 
          ncol = 1, nrow = 3, align = "hv", 
          widths = c(2, 1), heights = c(1, 2),
          common.legend = TRUE)


zplot <- ggdensity(data2$Middle, palette = "red")

ggarrange(xplot, zplot, 
          ncol = 1, nrow = 2,
          heights = c(1, 0.5))


ggarrange(xplot, yplot, zplot, 
          ncol = 1, nrow = 3,
          heights = c(1, 0.5, 0.3))




ggplot(histogram, aes(f0))  + 
  geom_histogram(data=data,fill = "red", alpha = 0.2) +
  geom_histogram(data=data2,fill = "blue", alpha = 0.2) # +
  # geom_histogram(data=subset(dat,yy == 'c'),fill = "green", alpha = 0.2)





data3 <- read.csv("breaks_comp.csv", sep=';')


plot_multi_histogram <- function(df, feature, label_column) {
  plt <- ggplot(df, aes(x=eval(parse(text=feature)), fill=eval(parse(text=label_column)))) +
    geom_histogram(alpha=0.7, position="identity", aes(y = ..density..), color="black") +
    geom_density(alpha=0.7) +
    geom_vline(aes(xintercept=mean(eval(parse(text=feature)))), color="black", linetype="dashed", size=1) +
    labs(x=feature, y = "Density")
  plt + guides(fill=guide_legend(title=label_column))
}
plot_multi_histogram(data3, 'Middle','Group')


