rm(list=ls(all=TRUE))
if (!require(tidyverse))    install.packages("tidyverse")
if (!require(furrr))        install.packages("furrr")
if (!require(here))         install.packages("here")
if (!require(raster))       install.packages("raster")
if (!require(cowplot))      install.packages("cowplot")
if (!require(skimr))        install.packages("skimr")
if (!require(pspearman))    install.packages("pspearman")
if (!require(ggstatsplot))  install.packages("ggstatsplot")
if (!require(ggasym))       install.packages("ggasym")


library(tidyverse)
library(here)
library(raster)
plots_dir <- here("Body/4Figures") %>% normalizePath() # Windows-people care

here("Body/4Figures/SlipAndJump.R.01.pdf") %>% normalizePath() %>% pdf()
  
##### 1: READ microhomology from pair-wise alignments
homol = read.table(here("Body/2Derived/HeatMaps/100x100.csv") %>% 
                     normalizePath(),
                   sep = ';',
                   header = TRUE)
row.names(homol) = homol$X
homol = homol[, -1]
# make long vertical table from the matrix
for (i in 1:nrow(homol)) {
  for (j in 1:ncol(homol)) {
    # i  = 2; j = 1
    FirstWindow = as.character(row.names(homol)[i])
    SecondWindow = as.character(names(homol)[j])
    Score = as.numeric(homol[i, j])
    OneLine = data.frame(FirstWindow, SecondWindow, Score)
    if (i == 1 & j == 1) {
      Final = OneLine
    }
    if (i > 1 | j > 1) {
      Final = rbind(Final, OneLine)
    }
  }
}

## the matrix is symmetric - I need to keep only one triangle: X>Y (don't need also diagonal, which is made by '500's)
Final$FirstWindow = as.numeric(as.character(Final$FirstWindow))
Final$SecondWindow = gsub('X', '', Final$SecondWindow) %>% as.numeric()

nrow(Final)
Final = Final[Final$FirstWindow > Final$SecondWindow, ]
nrow(Final)
MicroHomology = Final

##### 2: READ density of direct repeats per window

DirectRepDensity <- 
  here("Body/2Derived/HeatMaps/Link_matrix_direct_major_activ_left.csv") %>% 
  normalizePath() %>%
  read.table(sep = ';',
             header = TRUE)
DirectRepDensity = DirectRepDensity[, -1]

# make long vertical table from the matrix
for (i in 1:nrow(DirectRepDensity)) {
  for (j in 1:ncol(DirectRepDensity)) {
    # i  = 2; j = 1
    FirstWindow = as.character(row.names(DirectRepDensity)[i])
    SecondWindow = as.character(names(DirectRepDensity)[j])
    Score = as.numeric(DirectRepDensity[i, j])
    OneLine = data.frame(FirstWindow, SecondWindow, Score)
    if (i == 1 & j == 1) {
      Final = OneLine
    }
    if (i > 1 | j > 1) {
      Final = rbind(Final, OneLine)
    }
  }
}

## the matrix is symmetric - I need to keep only one triangle: X>Y (don't need also diagonal, which is made by '500's)
Final$FirstWindow = as.numeric(as.character(Final$FirstWindow))
Final$SecondWindow = gsub('X', '', Final$SecondWindow) %>% as.numeric()

nrow(Final)
Final = Final[Final$FirstWindow > Final$SecondWindow, ]
nrow(Final)
DirectRepDensity = Final

##### 3: correlate MicroHomology and  DirectRepDensity, derive HomologyAndRepeats dataset
DirectRepDensity = DirectRepDensity[
  order(DirectRepDensity$FirstWindow, DirectRepDensity$SecondWindow), 
  ]

MicroHomology = MicroHomology[
  order(MicroHomology$FirstWindow, MicroHomology$SecondWindow), 
  ]
pspearman::spearman.test(DirectRepDensity$Score, 
                         MicroHomology$Score) 
# S = 1.8841e+10, p-value = 1.707e-06 rho = 0.06796994 

nrow(DirectRepDensity) # 4950
plot(DirectRepDensity$Score, 
     MicroHomology$Score)

pspearman::spearman.test(DirectRepDensity[DirectRepDensity$Score > 0, ]$Score, 
                         MicroHomology[DirectRepDensity$Score > 0, ]$Score) 
# S = 1424127962, p-value = 0.006126 rho = 0.05998883 
plot(DirectRepDensity[DirectRepDensity$Score > 0, ]$Score, 
     MicroHomology[DirectRepDensity$Score > 0, ]$Score)

HomologyAndRepeats = data.frame(DirectRepDensity$FirstWindow,
                                DirectRepDensity$SecondWindow,
                                DirectRepDensity$Score,
                                MicroHomology$Score)
names(HomologyAndRepeats) = c('FirstWindow',
                              'SecondWindow',
                              'DirectRepeatsDensity',
                              'MicroHomologyScore')

##### 4: READ MITOBREAK AND FILTER (KEEP ONLY MAJOR ARC DELETIONS):

breaks = read.table(here("Body/1Raw/MitoBreakDB_12122019.csv") %>% 
                      normalizePath(), 
                    sep = ',', 
                    header = TRUE)
breaks$X5..breakpoint = as.numeric(as.character(breaks$X5..breakpoint))
summary(breaks$X5..breakpoint)
breaks$X3..breakpoint = as.numeric(as.character(breaks$X3..breakpoint))
summary(breaks$X3..breakpoint)
breaks = breaks[!is.na(breaks$X3..breakpoint) &
                  !is.na(breaks$X5..breakpoint), ]

breaks$FirstWindowBreakpoint = breaks$X3..breakpoint
breaks$SecondWindowBreakpoint = breaks$X5..breakpoint
breaks = breaks[breaks$FirstWindowBreakpoint   > 5781 &
                  breaks$FirstWindowBreakpoint < 16569 &
                  breaks$SecondWindowBreakpoint > 5781 &
                  breaks$SecondWindowBreakpoint < 16569, ] # can make it better!! take in ot account 0-100?
# поскольку координаты не такие простые (см ниже) - чтобы не париться можно взять все точки разрыва что больше чем 5781 и меньше чем 16569
# OH: 110-441
# OL: 5721-5781

HomologyAndRepeats$Deletion = 0
for (i in 1:nrow(HomologyAndRepeats)) {
  # i = 1
  FirstWindow  = HomologyAndRepeats$FirstWindow[i]
  SecondWindow = HomologyAndRepeats$SecondWindow[i]
  TempBreaks = breaks[breaks$FirstWindowBreakpoint   >= FirstWindow &
                        breaks$FirstWindowBreakpoint < (FirstWindow + 100) &
                        breaks$SecondWindowBreakpoint >= SecondWindow &
                        breaks$SecondWindowBreakpoint < (SecondWindow + 100), ]
  if (nrow(TempBreaks) > 0) {
    HomologyAndRepeats$Deletion[i] = 1
  }
}
  
table(HomologyAndRepeats$Deletion)
# 0    1 
# 4466  484 
skimr::skim(HomologyAndRepeats)

a <-
  glm(
    HomologyAndRepeats$Deletion ~ HomologyAndRepeats$DirectRepeatsDensity + HomologyAndRepeats$MicroHomologyScore,
    family = binomial
  )
summary(a)
ggstatsplot::ggcoefstats(a)

broom::tidy(a)
# term                            estimate std.error statistic  p.value
# <chr>                              <dbl>     <dbl>     <dbl>    <dbl>
# 1 (Intercept)                     -3.69      0.260     -14.2   1.12e-45
# 2 HomologyAndRepeats$DirectRepea…  0.00353   0.00605     0.583 5.60e- 1
# 3 HomologyAndRepeats$MicroHomolo…  0.0158    0.00273     5.78  7.63e- 9

# broom::augment(a)

broom::glance(a)
# null.deviance df.null logLik   AIC   BIC deviance df.residual
# <dbl>   <int>  <dbl> <dbl> <dbl>    <dbl>       <int>
# 3170.    4949 -1568. 3142. 3161.    3136.        4947



a <-
  glm(
    HomologyAndRepeats$Deletion ~ scale(HomologyAndRepeats$MicroHomologyScore),
    family = binomial
  )
summary(a) 
ggstatsplot::ggcoefstats(a)

broom::tidy(a)
# term                          estimate std.error statistic    p.value
# <chr>                            <dbl>     <dbl>     <dbl>      <dbl>
# 1 (Intercept)                     -2.25     0.0491    -45.8     0.     
# 2 scale(HomologyAndRepeats$Mic…    0.274    0.0470      5.84    5.13e-9

# broom::augment(a)
# 
broom::glance(a)
# null.deviance df.null logLik   AIC   BIC deviance df.residual
#         <dbl>   <int>  <dbl> <dbl> <dbl>    <dbl>       <int>
#   1     3170.    4949 -1568. 3140. 3153.    3136.        4948

nrow(HomologyAndRepeats)

### may be add perfect repeats of Orlov - yes or no for each deletion? and see, if microhomology still important?


##### 5: READ GLOBAL FOLDING:

GlobalFolding = here("Body/2Derived/HeatMaps/Link_matrix1000_major.csv") %>%
  normalizePath() %>%
  read.table(sep = ';',
             header = TRUE)
row.names(GlobalFolding) = GlobalFolding$X
GlobalFolding = GlobalFolding[, -1]

# make long vertical table from the matrix
for (i in 1:nrow(GlobalFolding)) {
  for (j in 1:ncol(GlobalFolding)) {
    # i  = 2; j = 1
    FirstWindow = as.character(row.names(GlobalFolding)[i])
    SecondWindow = as.character(names(GlobalFolding)[j])
    Score = as.numeric(GlobalFolding[i, j])
    OneLine = data.frame(FirstWindow, SecondWindow, Score)
    if (i == 1 & j == 1) {
      Final = OneLine
    }
    if (i > 1 | j > 1) {
      Final = rbind(Final, OneLine)
    }
  }
}

## the matrix is symmetric - I need to keep only one triangle: X>Y (don't need also diagonal, which is noizy and bold)
Final$FirstWindow = as.numeric(as.character(Final$FirstWindow))
Final$SecondWindow = gsub('X', '', Final$SecondWindow) %>% as.numeric()

nrow(Final)
Final = Final[Final$FirstWindow > Final$SecondWindow, ]
nrow(Final)
GlobalFolding1000 = Final
GlobalFolding1000 = GlobalFolding1000[order(GlobalFolding1000$FirstWindow,
                                            GlobalFolding1000$SecondWindow), ]
names(GlobalFolding1000) = c('FirstWindowWholeKbRes',
                             'SecondWindowWholeKbRes',
                             'GlobalFolding1000Score')

##### 5.1: READ GLOBAL FOLDING WITH WINDOW = 100 bp: (it automaticaly rewrites the GlobalFolding matrix from the previous point 5)

GlobalFolding = here("Body/2Derived/HeatMaps/Link_matrix100hydra_major.csv") %>%
  normalizePath() %>%
  read.table(sep = ';',
             header = TRUE)
row.names(GlobalFolding) = GlobalFolding$X
GlobalFolding = GlobalFolding[, -1]

# make long vertical table from the matrix
for (i in 1:nrow(GlobalFolding)) {
  for (j in 1:ncol(GlobalFolding)) {
    # i  = 2; j = 1
    FirstWindow = as.character(row.names(GlobalFolding)[i])
    SecondWindow = as.character(names(GlobalFolding)[j])
    Score = as.numeric(GlobalFolding[i, j])
    OneLine = data.frame(FirstWindow, SecondWindow, Score)
    if (i == 1 & j == 1) {
      Final = OneLine
    }
    if (i > 1 | j > 1) {
      Final = rbind(Final, OneLine)
    }
  }
}

## the matrix is symmetric - I need to keep only one triangle: X>Y (don't need also diagonal, which is made by '500's)
Final$FirstWindow = as.numeric(as.character(Final$FirstWindow))
Final$SecondWindow = gsub('X', '', Final$SecondWindow) %>% as.numeric()

nrow(Final)
Final = Final[Final$FirstWindow > Final$SecondWindow, ]
nrow(Final) 
# Should we delete bold diagonal or erase it to zeroes??? If delete, dimension will be decreased - try this. delete 5 windows next to diagonal (500)
nrow(Final)
Final = Final[Final$FirstWindow > Final$SecondWindow + 1000, ]
nrow(Final) # 500 or 1000!!!!!! similarly good results but 1000 is a bit better
GlobalFolding = Final
GlobalFolding = GlobalFolding[order(GlobalFolding$FirstWindow, GlobalFolding$SecondWindow), ]
names(GlobalFolding)[3] = c('GlobalFoldingScore')

pltHeatmap_gFolding_sw100 <- 
  here("Body/2Derived/HeatMaps/Link_matrix100hydra_major.csv") %>% 
  read.table(sep = ';', 
             header = TRUE) %>% 
  gather(-X, 
         key = "SecondWindow", 
         value = "Score") %>%
  rename(FirstWindow = X) %>% 
  mutate(SecondWindow = stringr::str_extract(SecondWindow, 
                                             "\\d+") %>% 
           as.integer()) %>%
  filter(FirstWindow  >= 5950,
         SecondWindow >= 5950) %>% 
  tibble() %>% 
  ggplot(aes(x = FirstWindow, 
             y = SecondWindow, 
             fill = Score)) +
  geom_tile() +
  scale_y_reverse() +
  scale_fill_viridis_c(option = "D", direction = -1) +
  theme_bw(base_size = 18)

cowplot::save_plot(
  plot = pltHeatmap_gFolding_sw100,
  base_height = 8.316,
  base_width = 11.594,
  file = normalizePath(
    file.path(plots_dir, 'heatmap_global_folding_sw100.pdf')
  )
)

# GlobalFolding - is the whole genome without bold diagonal, not only the major arc!! Keep only major arc in downstream analyses.
# will do it when merge with InvRepDens.

###### 6: READ INVERTED REPEATS WITH STEP 1000

InvRepDens = 
  here("Body/2Derived/HeatMaps/Link_matrix_1000_invert_major_activ_left.csv") %>%
  normalizePath() %>% 
  read.table(
  sep = ';',
  header = TRUE
)
row.names(InvRepDens) = InvRepDens$X
InvRepDens = InvRepDens[, -1]

# make long vertical table from the matrix
for (i in 1:nrow(InvRepDens)) {
  for (j in 1:ncol(InvRepDens)) {
    # i  = 2; j = 1
    FirstWindow = as.character(row.names(InvRepDens)[i])
    SecondWindow = as.character(names(InvRepDens)[j])
    Score = as.numeric(InvRepDens[i, j])
    OneLine = data.frame(FirstWindow, SecondWindow, Score)
    if (i == 1 & j == 1) {
      Final = OneLine
    }
    if (i > 1 | j > 1) {
      Final = rbind(Final, OneLine)
    }
  }
}

## the matrix is symmetric - I need to keep only one triangle: X>Y (don't need also diagonal, which is made by '500's)
Final$FirstWindow = as.numeric(as.character(Final$FirstWindow))
Final$SecondWindow = gsub('X', '', Final$SecondWindow) %>% as.numeric()

nrow(Final)
Final = Final[Final$FirstWindow > Final$SecondWindow,]
nrow(Final)
InvRepDens = Final
InvRepDens = InvRepDens[order(InvRepDens$FirstWindow, InvRepDens$SecondWindow),]

###### 6.1: READ INVERTED REPEATS WITH STEP 100 bp (it automatically rewrites InvRepDens from previous point 6)

InvRepDens = 
  here("Body/2Derived/HeatMaps/Link_matrix_invert_major_activ_left.modified.csv") %>% 
  normalizePath() %>% 
  read.table(
    sep = '\t',
  header = TRUE,
  row.names = 1
) # , row.names = NULL)

# make long vertical table from the matrix
for (i in 1:nrow(InvRepDens)) {
  for (j in 1:ncol(InvRepDens)) {
    # i  = 2; j = 1
    FirstWindow = as.character(row.names(InvRepDens)[i])
    SecondWindow = as.character(names(InvRepDens)[j])
    Score = as.numeric(InvRepDens[i, j])
    OneLine = data.frame(FirstWindow, SecondWindow, Score)
    if (i == 1 & j == 1) {
      Final = OneLine
    }
    if (i > 1 | j > 1) {
      Final = rbind(Final, OneLine)
    }
  }
}

## the matrix is symmetric - I need to keep only one triangle: X>Y (don't need also diagonal, which is made by '500's)
Final$FirstWindow = as.numeric(as.character(Final$FirstWindow))
Final$SecondWindow = gsub('X', '', Final$SecondWindow) %>% as.numeric()

nrow(Final)
Final = Final[Final$FirstWindow > Final$SecondWindow,]
nrow(Final)
InvRepDens = Final
InvRepDens = InvRepDens[order(InvRepDens$FirstWindow, InvRepDens$SecondWindow),]
summary(InvRepDens$FirstWindow)  # 6000 15800
summary(InvRepDens$SecondWindow) # 5900 15700
names(InvRepDens)[3] = c('InvRepDensScore')
skimr::skim(InvRepDens)

###### 6.2: READ INVERTED REPEATS WITH STEP 100 bp WITH OVERLAPS (it automatically rewrites InvRepDens from previous point 6.1)
# not nice results - use 6.1 
#InvRepDens = read.table("../../Body/2Derived/HeatMaps/Link_matrix_direct_major_activ_left_cross.KpModifByHand.mtrx", sep = '\t',header = TRUE, row.names = 1) # , row.names = NULL)

# make long vertical table from the matrix
#for (i in 1:nrow(InvRepDens))
#{
#  for (j in 1:ncol(InvRepDens))
#  { # i  = 2; j = 1
#    FirstWindow = as.character(row.names(InvRepDens)[i])
#    SecondWindow = as.character(names(InvRepDens)[j])
#    Score = as.numeric(InvRepDens[i,j])
#    OneLine = data.frame(FirstWindow,SecondWindow,Score)
#    if (i == 1 & j == 1) {Final = OneLine}
#    if (i > 1 | j > 1) {Final = rbind(Final,OneLine)}
#  }
#}

## the matrix is symmetric - I need to keep only one triangle: X>Y (don't need also diagonal)
#Final$SecondWindow = gsub('X','',Final$SecondWindow)
#Final$FirstWindow = as.numeric(as.character(Final$FirstWindow)); Final$SecondWindow = as.numeric(Final$SecondWindow); 
#nrow(Final); Final=Final[Final$FirstWindow > Final$SecondWindow,]; nrow(Final)  
#InvRepDens = Final
#InvRepDens = InvRepDens[order(InvRepDens$FirstWindow,InvRepDens$SecondWindow),]
#summary(InvRepDens$FirstWindow)  # 6000 15800
#summary(InvRepDens$SecondWindow) # 5900 15700
#names(InvRepDens)[3] = c('InvRepDensScore');

###### 7: CORRELATE GlobalFolding$Score and InvRepDens$Score - weak positive!
merged = merge(InvRepDens, GlobalFolding, by = c("FirstWindow", "SecondWindow"))
summary(merged$FirstWindow)  # diag 500: 6500 15800; diag 1000: 7000 - 15800
summary(merged$SecondWindow) # diag 500: 5900 15200; diag 1000: 5900 14700
skimr::skim(merged)
pspearman::spearman.test(merged$InvRepDensScore, merged$GlobalFoldingScore) 
# S = 1.0177e+10, p-value = 0.001745 rho = 0.04945796
# diag 500: rho = 0.04926082, p-value = 0.0009922; diag 1000: rho = 0.04945796, p = 0.001743
nrow(merged) # 4005

###### 8: ADD InfinitySign parameter into HomologyAndRepeats dataset (13 - 16 kb vs 6-9 kb):
HomologyAndRepeats$InfinitySign = 0
for (i in 1:nrow(HomologyAndRepeats)) {
  if (HomologyAndRepeats$FirstWindow[i] >= 13000 &
      HomologyAndRepeats$FirstWindow[i] <= 16000 &
      HomologyAndRepeats$SecondWindow[i] >= 6000 &
      HomologyAndRepeats$SecondWindow[i] <= 9000) {
    HomologyAndRepeats$InfinitySign[i] = 1
  }
}
janitor::tabyl(HomologyAndRepeats$InfinitySign)
janitor::tabyl(HomologyAndRepeats, Deletion, InfinitySign)
sjmisc::frq(HomologyAndRepeats)
HomologyAndRepeats %>% skim()
HomologyAndRepeats %>% group_by(Deletion) %>% skim()
HomologyAndRepeats %>% group_by(InfinitySign) %>% skim()

## merge HomologyAndRepeats with merged(InvRepDens + GlobalFolding)
dim(HomologyAndRepeats) # 4950
HomologyAndRepeats = merge(HomologyAndRepeats,
                           merged,
                           by = c("FirstWindow", "SecondWindow"))
dim(HomologyAndRepeats) # diag 500: 4465; diag 1000:  4005

# is GlobalFoldingScore higher within the cross according to our InfinitySign model? YES!!! 
wilcox.test(HomologyAndRepeats[HomologyAndRepeats$InfinitySign == 1, ]$GlobalFoldingScore,
            HomologyAndRepeats[HomologyAndRepeats$InfinitySign == 0, ]$GlobalFoldingScore)# diag 1000: 3.358e-09
boxplot(
  HomologyAndRepeats[HomologyAndRepeats$InfinitySign == 1, ]$GlobalFoldingScore,
  HomologyAndRepeats[HomologyAndRepeats$InfinitySign == 0, ]$GlobalFoldingScore,
  notch = TRUE,
  names = c('stem', 'loop'),
  ylab = 'in silico folding score',
  outline = FALSE
)

pltViolRepFoldingInfSign <- 
  ggbetweenstats(
    data = HomologyAndRepeats,
    x = InfinitySign,
    y = GlobalFoldingScore,
    type = "np", 
    # Wilcoxon for two group
    mean.ci = TRUE, 
    nboot = 10000,
    # number of iteration for statistical CI
    k = 5,
    # number of decimal places for statistical results
    outlier.tagging = TRUE, # whether outliers need to be tagged
    outlier.label = Deletion,
    xlab = '"3D" Position',
    # label for the x-axis variable
    ylab = "in silico folding score",
    # label for the y-axis variable
    title = "The effect of repeats' position on folding score",
    # title text for the plot
    ggtheme = ggthemes::theme_fivethirtyeight(),
    # choosing a different theme
    package = "wesanderson",
    # package from which color palette is to be taken
    palette = "Royal1",
    # choosing a different color palette
    notch = TRUE,
    messages = TRUE
  )
# Note: 95% CI for effect size estimate was computed with 10000 bootstrap samples
# Note: Shapiro-Wilk Normality Test for in silico folding score: p-value = < 0.001
# Note: Bartlett's test for homogeneity of variances for factor "3D" Position: p-value = < 0.001
cowplot::save_plot(
  plot = pltViolRepFoldingInfSign,
  base_height = 8,
  base_asp = 1.618,
  file = normalizePath(
    file.path(plots_dir, 'violin_rep_folding_infsign_np.pdf')
  )
)



t.test(HomologyAndRepeats[HomologyAndRepeats$InfinitySign == 1, ]$GlobalFoldingScore,
       HomologyAndRepeats[HomologyAndRepeats$InfinitySign == 0, ]$GlobalFoldingScore) # diag 1000: 0.002639

pltViolRepFoldingInfSignP <- 
  ggbetweenstats(
    data = HomologyAndRepeats,
    x = InfinitySign,
    y = GlobalFoldingScore,
    type = "p", 
    # Wilcoxon for two group
    mean.ci = TRUE, 
    nboot = 10000,
    # number of iteration for statistical CI
    k = 5,
    # number of decimal places for statistical results
    outlier.tagging = TRUE, # whether outliers need to be tagged
    outlier.label = FirstWindow,
    xlab = '"3D" Position',
    # label for the x-axis variable
    ylab = "in silico folding score",
    # label for the y-axis variable
    title = "The effect of repeats' position on folding score",
    # title text for the plot
    ggtheme = ggthemes::theme_fivethirtyeight(),
    # choosing a different theme
    package = "wesanderson",
    # package from which color palette is to be taken
    palette = "Royal1",
    # choosing a different color palette
    notch = TRUE,
    messages = TRUE
  )
# Note: 95% CI for effect size estimate was computed with 10000 bootstrap samples
# Note: Shapiro-Wilk Normality Test for in silico folding score: p-value = < 0.001
# Note: Bartlett's test for homogeneity of variances for factor "3D" Position: p-value = < 0.001
cowplot::save_plot(
  plot = pltViolRepFoldingInfSignP,
  base_height = 8,
  base_asp = 1.618,
  file = normalizePath(
    file.path(plots_dir, 'violin_rep_folding_infsign_p.pdf')
  )
)


summary(HomologyAndRepeats[HomologyAndRepeats$InfinitySign == 1, ]$GlobalFoldingScore) # diag 1000: 0.3893
summary(HomologyAndRepeats[HomologyAndRepeats$InfinitySign == 0, ]$GlobalFoldingScore) # diag 1000: 0.09755

# we have to link better global folding and InfinitySign model - till now it was done by eye. Clusterisation? One cluster? 
# dev.off()
###### 9: LOGISTIC REGRESSION: HomologyAndRepeats$Deletion as a function of HomologyAndRepeats$MicroHomologyScore and HomologyAndRepeats$InfinitySign:

a <-
  glm(
    HomologyAndRepeats$Deletion ~ HomologyAndRepeats$MicroHomologyScore + HomologyAndRepeats$InfinitySign,
    family = 'binomial'
  )
summary(a)
ggstatsplot::ggcoefstats(a)
broom::tidy(a)
# term                              estimate std.error statistic  p.value
# <chr>                                <dbl>     <dbl>     <dbl>    <dbl>
#   1 (Intercept)                      -4.57     0.298      -15.4  2.52e-53
# 2 HomologyAndRepeats$MicroHomology…   0.0189   0.00303     6.26  3.80e-10
# 3 HomologyAndRepeats$InfinitySign     2.17     0.106       20.4  1.24e-92
# broom::augment(a)
broom::glance(a)
# null.deviance df.null logLik   AIC   BIC deviance df.residual
#          <dbl>   <int>  <dbl> <dbl> <dbl>    <dbl>       <int>
#          2925.    4004 -1225. 2457. 2475.    2451.        4002

a <-
  glm(
    HomologyAndRepeats$Deletion ~ scale(HomologyAndRepeats$MicroHomologyScore) + scale(HomologyAndRepeats$InfinitySign),
    family = 'binomial'
  )
summary(a) # PAPER!!! 0.33 + 0.91
ggstatsplot::ggcoefstats(a)
broom::tidy(a)
# term                             estimate std.error statistic   p.value
# <chr>                               <dbl>     <dbl>     <dbl>     <dbl>
# 1 (Intercept)                        -2.38     0.0641    -37.2  2.37e-302
# 2 scale(HomologyAndRepeats$MicroH…    0.326    0.0520     6.26  3.80e- 10
# 3 scale(HomologyAndRepeats$Infini…    0.906    0.0444     20.4  1.24e- 92

# broom::augment(a)
broom::glance(a)
# null.deviance df.null logLik   AIC   BIC deviance df.residual
#          <dbl>   <int>  <dbl> <dbl> <dbl>    <dbl>       <int>
#          2925.    4004 -1225. 2457. 2475.    2451.        4002

a <-
  glm(
    HomologyAndRepeats$Deletion ~ HomologyAndRepeats$MicroHomologyScore + HomologyAndRepeats$GlobalFoldingScore,
    family = 'binomial'
  )
summary(a) # non significant - may be I have to take it on bigger scale! (1kb without diagonal, because this is global parameter not precise)
# to reconstruct 100 bp matrix back from 1 kb matrix!!!!!
ggstatsplot::ggcoefstats(a)
broom::tidy(a)
# term                              estimate std.error statistic  p.value
# <chr>                                <dbl>     <dbl>     <dbl>    <dbl>
# 1 (Intercept)                       -3.57      0.267     -13.4   6.91e-41
# 2 HomologyAndRepeats$MicroHomology…  0.0171    0.00279     6.12  9.54e-10
# 3 HomologyAndRepeats$GlobalFolding…  0.00558   0.0288      0.194 8.46e- 1
# broom::augment(a)
broom::glance(a)
# null.deviance df.null logLik   AIC   BIC deviance df.residual
#          <dbl>   <int>  <dbl> <dbl> <dbl>    <dbl>       <int>
#          2925.    4004 -1444.2893. 2912.    2887.        4002

# get residuals and correlate them with global matrix

a <-
  glm(HomologyAndRepeats$Deletion ~ HomologyAndRepeats$MicroHomologyScore,
      family = 'binomial')
ggstatsplot::ggcoefstats(a)
# broom::augment(a) %>% dplyr::select("HomologyAndRepeats.Deletion", "HomologyAndRepeats.MicroHomologyScore", ".resid") %>% ggplot(aes(x = .resid, y = HomologyAndRepeats.MicroHomologyScore, colour = HomologyAndRepeats.Deletion), data = .) + geom_point()
HomologyAndRepeats$Residuals = residuals(a)
summary(HomologyAndRepeats$Residuals)

pspearman::spearman.test(HomologyAndRepeats$Residuals,
                         HomologyAndRepeats$GlobalFoldingScore) # S = 1.0104e+10, p-value = 0.0003666 rho = 0.05628272 
# diag 500: rho = 0.04714374, p = 0.001627; diag 1000: rho = 0.05628272, p = 0.0003658

#### reconstruct Global folding 100 bp back from 1kb resolution (GlobalFolding1000) assuming that global folding can work remotely enough.
#### another idea - to use a distance from a given cell to closest contact (from global matrix) - so, infinity sign is not zero or one, but continuos!

# round(6600,-3) = 7000; round(6500,-3) = 6000;
HomologyAndRepeats$FirstWindowWholeKbRes = round(HomologyAndRepeats$FirstWindow, -3)
HomologyAndRepeats$SecondWindowWholeKbRes = round(HomologyAndRepeats$SecondWindow, -3)
table(HomologyAndRepeats$FirstWindowWholeKbRes)
table(HomologyAndRepeats$SecondWindowWholeKbRes)
janitor::tabyl(HomologyAndRepeats, 
               FirstWindowWholeKbRes, 
               SecondWindowWholeKbRes)

nrow(HomologyAndRepeats)  # diag 1000: 4005
HomologyAndRepeats = merge(
  HomologyAndRepeats,
  GlobalFolding1000,
  by = c("FirstWindowWholeKbRes", "SecondWindowWholeKbRes")
)
nrow(HomologyAndRepeats)  # diag 1000: 4005

a <-
  glm(
    HomologyAndRepeats$Deletion ~ scale(HomologyAndRepeats$MicroHomologyScore) + scale(HomologyAndRepeats$GlobalFolding1000Score),
    family = 'binomial'
  )
summary(a)
ggstatsplot::ggcoefstats(a)
broom::tidy(a) %>% knitr::kable(format = "markdown")
# |term                                             |   estimate| std.error|  statistic|   p.value|
# |:------------------------------------------------|----------:|---------:|----------:|---------:|
# |(Intercept)                                      | -2.0372126| 0.0503135| -40.490337| 0.0000000|
# |scale(HomologyAndRepeats$MicroHomologyScore)     |  0.2912598| 0.0480566|   6.060761| 0.0000000|
# |scale(HomologyAndRepeats$GlobalFolding1000Score) |  0.0912813| 0.0444382|   2.054118| 0.0399643|
# broom::augment(a)
broom::glance(a) %>% knitr::kable(format = "markdown")
# | null.deviance| df.null|    logLik|      AIC|      BIC| deviance| df.residual|
# |-------------:|-------:|---------:|--------:|--------:|--------:|-----------:|
# |      2924.693|    4004| -1441.768| 2889.536| 2908.422| 2883.536|        4002|

# diag 500:
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                                      -2.13660    0.05004 -42.698  < 2e-16 ***
# scale(HomologyAndRepeats$MicroHomologyScore)      0.29139    0.04774   6.103 1.04e-09 ***
# scale(HomologyAndRepeats$GlobalFolding1000Score)  0.07605    0.04477   1.699   0.0894 .  


##### derive distance to the strongest contact: 6500 vs 14500 (see heatmap: global folding 1 kb resolution). Check coordinates once more!!!
HomologyAndRepeats$DistanceToContact = 0
for (i in 1:nrow(HomologyAndRepeats)) {
  # i = 1
  HomologyAndRepeats$DistanceToContact[i] = raster::pointDistance(
    c(
      HomologyAndRepeats$FirstWindow[i],
      HomologyAndRepeats$SecondWindow[i]
    ),
    c(14550, 6550),
    lonlat = FALSE
  )
}

skim(HomologyAndRepeats$DistanceToContact) %>% knitr::kable()
# |skim_type |skim_variable | n_missing| complete_rate| numeric.mean| numeric.sd| numeric.p0| numeric.p25| numeric.p50| numeric.p75| numeric.p100|numeric.hist |
# |:---------|:-------------|---------:|-------------:|------------:|----------:|----------:|-----------:|-----------:|-----------:|------------:|:------------|
# |numeric   |data          |         0|             1|     3726.984|   1747.957|   70.71068|    2392.697|    3842.525|    4962.358|     8245.302|▃▆▇▅▁        |
summary(HomologyAndRepeats$DistanceToContact) # the closest: -70; the most distant: -8245


a <-
  glm(
    HomologyAndRepeats$Deletion ~ scale(HomologyAndRepeats$MicroHomologyScore) + scale(HomologyAndRepeats$DistanceToContact),
    family = 'binomial'
  )
summary(a)
ggstatsplot::ggcoefstats(a) # Error in grid.Call(C_convert, x, as.integer(whatfrom), as.integer(whatto),  : 
# Viewport has zero dimension(s) !!!!
broom::tidy(a) %>% knitr::kable()
# |term                                         |  estimate| std.error| statistic| p.value|
# |:--------------------------------------------|---------:|---------:|---------:|-------:|
# |(Intercept)                                  | -2.566536| 0.0738349| -34.76046|       0|
# |scale(HomologyAndRepeats$MicroHomologyScore) |  0.432243| 0.0526914|   8.20330|       0|
# |scale(HomologyAndRepeats$DistanceToContact)  | -1.248810| 0.0652043| -19.15226|       0|
# broom::augment(a)
broom::glance(a) %>% knitr::kable()
# | null.deviance| df.null|    logLik|      AIC|      BIC| deviance| df.residual|
# |-------------:|-------:|---------:|--------:|--------:|--------:|-----------:|
# |      2924.693|    4004| -1201.136| 2408.272| 2427.158| 2402.272|        4002|

#Coefficients:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)                                    -2.56654    0.07383 -34.760  < 2e-16 ***
#  scale(HomologyAndRepeats$MicroHomologyScore)  0.43224    0.05269   8.203 2.34e-16 ***
#  scale(HomologyAndRepeats$DistanceToContact)  -1.24881    0.06520 -19.152  < 2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Dispersion parameter for binomial family taken to be 1)
# Null deviance: 2924.7  on 4004  degrees of freedom
# Residual deviance: 2402.3  on 4002  degrees of freedom
# AIC: 2408.3

### derive distance to common repeat:
HomologyAndRepeats$DistanceToContact = 0
for (i in 1:nrow(HomologyAndRepeats)) {
  # i = 1
  HomologyAndRepeats$DistanceToContact[i] = pointDistance(
    c(
      HomologyAndRepeats$FirstWindow[i],
      HomologyAndRepeats$SecondWindow[i]
    ),
    c(13447, 8469),
    lonlat = FALSE
  )  #  (8469-8482 - 13447-13459)
}
skim(HomologyAndRepeats$DistanceToContact) %>% knitr::kable()
# |skim_type |skim_variable | n_missing| complete_rate| numeric.mean| numeric.sd| numeric.p0| numeric.p25| numeric.p50| numeric.p75| numeric.p100|numeric.hist |
# |:---------|:-------------|---------:|-------------:|------------:|----------:|----------:|-----------:|-----------:|-----------:|------------:|:------------|
# |numeric   |data          |         0|             1|     2744.767|   1354.259|   56.30275|    1783.359|     2535.62|    3547.897|     6939.998|▃▇▅▂▁        |
summary(HomologyAndRepeats$DistanceToContact) # the closest: -56.3; the most distant: -6940.0 
a <-
  glm(
    HomologyAndRepeats$Deletion ~ scale(HomologyAndRepeats$MicroHomologyScore) + scale(HomologyAndRepeats$DistanceToContact),
    family = 'binomial'
  )
summary(a)
ggstatsplot::ggcoefstats(a) # Error in grid.Call(C_convert, x, as.integer(whatfrom), as.integer(whatto),  : 
# Viewport has zero dimension(s)
broom::tidy(a) %>% knitr::kable()
# |term                                         |   estimate| std.error|  statistic| p.value|
# |:--------------------------------------------|----------:|---------:|----------:|-------:|
# |(Intercept)                                  | -2.3579307| 0.0653587| -36.076742|       0|
# |scale(HomologyAndRepeats$MicroHomologyScore) |  0.2766745| 0.0503802|   5.491729|       0|
# |scale(HomologyAndRepeats$DistanceToContact)  | -1.0060227| 0.0687794| -14.626799|       0|
# broom::augment(a)
broom::glance(a) %>% knitr::kable()
# | null.deviance| df.null|    logLik|      AIC|     BIC| deviance| df.residual|
# |-------------:|-------:|---------:|--------:|-------:|--------:|-----------:|
# |      2924.693|    4004| -1304.572| 2615.144| 2634.03| 2609.144|        4002|

#                                             Estimate Std. Error z value Pr(>|z|)    
#(Intercept)                                  -2.35793    0.06536 -36.077  < 2e-16 ***
#  scale(HomologyAndRepeats$MicroHomologyScore)  0.27667    0.05038   5.492 3.98e-08 ***
#  scale(HomologyAndRepeats$DistanceToContact)  -1.00602    0.06878 -14.627  < 2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for binomial family taken to be 1)
#
#Null deviance: 2924.7  on 4004  degrees of freedom
#Residual deviance: 2609.1  on 4002  degrees of freedom
#AIC: 2615.1

dev.off()  


###### run many log regr in the loop and find the Contact Point with the best AIC 
library(furrr) # for parallel quantification
plan("multiprocess")

mod_fun <- function(Coord1, Coord2) {
  TempDF <- 
    HomologyAndRepeats %>% 
    dplyr::select(FirstWindow, 
           SecondWindow, 
           Deletion, 
           MicroHomologyScore) %>% 
    mutate(TempDistanceToContact = purrr::map2_dbl(FirstWindow, 
                                               SecondWindow, 
                                               ~raster::pointDistance(c(.x, .y),
                                                                      c(Coord1,
                                                                        Coord2),
                                                                      lonlat = F)
                                               )
    )
  a <-
    glm(
      TempDF$Deletion ~ scale(TempDF$MicroHomologyScore) + scale(TempDF$TempDistanceToContact),
      family = 'binomial'
    )
  return(a)
}

Results <- HomologyAndRepeats %>%
  mutate(
    LogRegr.ContactPoint.Coord1 = FirstWindow + 50,
    LogRegr.ContactPoint.Coord2 = SecondWindow + 50
  ) %>%
  mutate(model = furrr::future_map2(LogRegr.ContactPoint.Coord1, 
                                    LogRegr.ContactPoint.Coord2, 
                                    mod_fun))

ResultsTidy <- Results$model %>% purrr::map_dfr( ~ broom::tidy(.x)[3,])
ResultsGlance <- Results$model %>% purrr::map_dfr(broom::glance)

HomologyAndRepeats$LogRegr.ContactPoint.PiValue <- ResultsTidy$p.value
HomologyAndRepeats$LogRegr.ContactPoint.Coeff <- ResultsTidy$estimate
HomologyAndRepeats$LogRegr.ContactPoint.AIC <- ResultsGlance$AIC
HomologyAndRepeats$LogRegr.ContactPoint.ResidualDeviance <- ResultsGlance$deviance
HomologyAndRepeats$LogRegr.ContactPoint.Coord1 <- Results$LogRegr.ContactPoint.Coord1
HomologyAndRepeats$LogRegr.ContactPoint.Coord2 <- Results$LogRegr.ContactPoint.Coord2

rm(Results, ResultsGlance, ResultsTidy)
invisible(gc())

write.table(
  HomologyAndRepeats,
  here("Body/3Results/SlipAndJump.HomologyAndRepeats.txt") %>% 
    normalizePath(),
  sep = '\t'
)

HomologyAndRepeats = read.table(
  here("Body/3Results/SlipAndJump.HomologyAndRepeats.txt") %>% 
    normalizePath(),
  sep = '\t'
  )

HomologyAndRepeats = HomologyAndRepeats[order(HomologyAndRepeats$LogRegr.ContactPoint.AIC), ]
names(HomologyAndRepeats)
summary(HomologyAndRepeats$LogRegr.ContactPoint.ResidualDeviance)
summary(HomologyAndRepeats$LogRegr.ContactPoint.AIC)
HomologyAndRepeats %>% skim()

temp = HomologyAndRepeats[
  HomologyAndRepeats$LogRegr.ContactPoint.Coord1 == 11950 & 
    HomologyAndRepeats$LogRegr.ContactPoint.Coord2 == 8950, 
  ]
temp


here("Body/4Figures/SlipAndJump.R.02.pdf") %>% normalizePath() %>% pdf()

par(mfrow = c(2, 4))

plot(
  HomologyAndRepeats$LogRegr.ContactPoint.Coord2,
  HomologyAndRepeats$LogRegr.ContactPoint.AIC,
  xlab = '5 prime position',
  ylab = 'AIC'
)
abline(v = 9000, col = 'red', lwd = 1)
abline(v = 6000, col = 'red', lwd = 1)
ggscatterstats(
  data = HomologyAndRepeats,
  x = LogRegr.ContactPoint.Coord2,
  y = LogRegr.ContactPoint.AIC,
  type = "robust",
  conf.level = 0.95,
  nboot = 10000,
  k = 5,
  marginal.type = "density",
  ggtheme = ggthemes::theme_fivethirtyeight(),
  # choosing a different theme
  messages = TRUE
)


plot(
  HomologyAndRepeats$LogRegr.ContactPoint.Coord2,
  HomologyAndRepeats$LogRegr.ContactPoint.Coeff,
  xlab = '5 prime position',
  ylab = 'Coefficient'
)
abline(v = 9000, col = 'red', lwd = 1)
abline(v = 6000, col = 'red', lwd = 1)
ggscatterstats(
  data = HomologyAndRepeats,
  x = LogRegr.ContactPoint.Coord2,
  y = LogRegr.ContactPoint.Coeff,
  type = "robust",
  conf.level = 0.95,
  nboot = 10000,
  k = 5,
  marginal.type = "density",
  ggtheme = ggthemes::theme_fivethirtyeight(),
  # choosing a different theme
  messages = TRUE
)


plot(
  HomologyAndRepeats$LogRegr.ContactPoint.Coord2,
  -log10(HomologyAndRepeats$LogRegr.ContactPoint.PiValue),
  xlab = '5 prime position',
  ylab = '-log10(p-value)'
)
abline(v = 9000, col = 'red', lwd = 1)
abline(v = 6000, col = 'red', lwd = 1)
ggscatterstats(
  data = HomologyAndRepeats,
  x = LogRegr.ContactPoint.Coord2,
  y = LogRegr.ContactPoint.PiValue,
  type = "robust",
  conf.level = 0.95,
  nboot = 10000,
  k = 5,
  marginal.type = "density",
  ggtheme = ggthemes::theme_fivethirtyeight(),
  # choosing a different theme
  messages = TRUE
)


plot(
  HomologyAndRepeats$LogRegr.ContactPoint.Coord2,
  HomologyAndRepeats$LogRegr.ContactPoint.ResidualDeviance,
  xlab = '5 prime position',
  ylab = 'ResidualDeviance'
)
abline(v = 9000, col = 'red', lwd = 1)
abline(v = 6000, col = 'red', lwd = 1) 
ggscatterstats(
  data = HomologyAndRepeats,
  x = LogRegr.ContactPoint.Coord2,
  y = LogRegr.ContactPoint.ResidualDeviance,
  type = "robust",
  conf.level = 0.95,
  nboot = 10000,
  k = 5,
  marginal.type = "density",
  ggtheme = ggthemes::theme_fivethirtyeight(),
  # choosing a different theme
  messages = TRUE
)

plot(
  HomologyAndRepeats$LogRegr.ContactPoint.Coord1,
  HomologyAndRepeats$LogRegr.ContactPoint.AIC,
  xlab = '3 prime position',
  ylab = 'AIC'
)
abline(v = 13000, col = 'red', lwd = 1)
abline(v = 16000, col = 'red', lwd = 1)
ggscatterstats(
  data = HomologyAndRepeats,
  x = LogRegr.ContactPoint.Coord1,
  y = LogRegr.ContactPoint.AIC,
  type = "robust",
  conf.level = 0.95,
  nboot = 10000,
  k = 5,
  marginal.type = "density",
  ggtheme = ggthemes::theme_fivethirtyeight(),
  # choosing a different theme
  messages = TRUE
)

plot(
  HomologyAndRepeats$LogRegr.ContactPoint.Coord1,
  HomologyAndRepeats$LogRegr.ContactPoint.Coeff,
  xlab = '3 prime position',
  ylab = 'Coefficient'
)
abline(v = 13000, col = 'red', lwd = 1)
abline(v = 16000, col = 'red', lwd = 1)
ggscatterstats(
  data = HomologyAndRepeats,
  x = LogRegr.ContactPoint.Coord1,
  y = LogRegr.ContactPoint.Coeff,
  type = "robust",
  conf.level = 0.95,
  nboot = 10000,
  k = 5,
  marginal.type = "density",
  ggtheme = ggthemes::theme_fivethirtyeight(),
  # choosing a different theme
  messages = TRUE
)


plot(
  HomologyAndRepeats$LogRegr.ContactPoint.Coord1,
  -log10(HomologyAndRepeats$LogRegr.ContactPoint.PiValue),
  xlab = '3 prime position',
  ylab = '-log10(p-value)'
)
abline(v = 13000, col = 'red', lwd = 1)
abline(v = 16000, col = 'red', lwd = 1)
ggscatterstats(
  data = HomologyAndRepeats,
  x = LogRegr.ContactPoint.Coord1,
  y = LogRegr.ContactPoint.PiValue,
  type = "robust",
  conf.level = 0.95,
  nboot = 10000,
  k = 5,
  marginal.type = "density",
  ggtheme = ggthemes::theme_fivethirtyeight(),
  # choosing a different theme
  messages = TRUE
)


plot(
  HomologyAndRepeats$LogRegr.ContactPoint.Coord1,
  HomologyAndRepeats$LogRegr.ContactPoint.ResidualDeviance,
  xlab = '3 prime position',
  ylab = 'ResidualDeviance'
)
abline(v = 13000, col = 'red', lwd = 1)
abline(v = 16000, col = 'red', lwd = 1) 
ggscatterstats(
  data = HomologyAndRepeats,
  x = LogRegr.ContactPoint.Coord1,
  y = LogRegr.ContactPoint.ResidualDeviance,
  type = "robust",
  conf.level = 0.95,
  nboot = 10000,
  k = 5,
  marginal.type = "density",
  ggtheme = ggthemes::theme_fivethirtyeight(),
  # choosing a different theme
  messages = TRUE
)


dev.off()  

##### Heatmap merged microhomology and AIC scores #### 
tib <- HomologyAndRepeats %>% 
  dplyr::select(
    LogRegr.ContactPoint.Coord1,
    LogRegr.ContactPoint.Coord2,
    MicroHomologyScore,
    LogRegr.ContactPoint.AIC
  ) %>% 
  ggasym::asymmetrise(., 
              LogRegr.ContactPoint.Coord1, 
              LogRegr.ContactPoint.Coord2)

pltHeatmap_mhAIC <- ggplot(tib,
                           aes(x = LogRegr.ContactPoint.Coord1, 
                               y = LogRegr.ContactPoint.Coord2)) +
  geom_asymmat(aes(fill_br = LogRegr.ContactPoint.AIC,
                   fill_tl = MicroHomologyScore)) +
  scale_fill_br_distiller(
    type = "seq",
    palette = "Spectral",
    direction = 1,
    na.value = "white",
    guide = guide_colourbar(
      direction = "horizontal",
      order = 1,
      title.position = "top"
    )
  ) +
  scale_fill_tl_distiller(
    type = "seq",
    palette = "RdYlGn",
    direction = 1,
    na.value = "white",
    guide = guide_colourbar(
      direction = "horizontal",
      order = 2,
      title.position = "top"
    )
  ) +
  labs(fill_br = "top-right Contact AIC",
       fill_tl = "bottom-left Microhomology",
       title = "Model of mtDNA contacts") +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank()
  )

cowplot::save_plot(
  plot = pltHeatmap_mhAIC,
  base_height = 8.316,
  base_width = 11.594,
  file = normalizePath(
    file.path(plots_dir, 'heatmap_microhomology_AIC.pdf')
  )
)

##### Again: Heatmap merged microhomology and AIC scores with actual deledions circles ####
tib <- HomologyAndRepeats %>% 
  dplyr::select(
    LogRegr.ContactPoint.Coord1,
    LogRegr.ContactPoint.Coord2,
    MicroHomologyScore,
    LogRegr.ContactPoint.AIC,
    Deletion
  ) %>% 
  asymmetrise(., 
              LogRegr.ContactPoint.Coord1, 
              LogRegr.ContactPoint.Coord2)

pltHeatmap_mhAIC_wt_deletions <- ggplot(tib,
                                        aes(x = LogRegr.ContactPoint.Coord1, 
                                            y = LogRegr.ContactPoint.Coord2)) +
  geom_asymmat(aes(fill_br = LogRegr.ContactPoint.AIC,
                   fill_tl = MicroHomologyScore)) +
  geom_point(
    data = subset(tib, Deletion > 0),
    aes(alpha = 0.3),
    shape = 1,
    size = 2.3,
    color = "#041c00",
    show.legend = FALSE,
    na.rm = TRUE
  ) +
  scale_fill_br_gradientn(
    colours = viridis::viridis(200),
    na.value = "white",
    guide = guide_colourbar(
      direction = "horizontal",
      order = 1,
      title.position = "top"
    )
  ) +
  scale_fill_tl_distiller(
    type = "seq",
    palette = "RdYlGn",
    direction = 1,
    na.value = "white",
    guide = guide_colourbar(
      direction = "horizontal",
      order = 2,
      title.position = "top"
    )
  ) +
  labs(fill_br = "top-right Contact AIC",
       fill_tl = "bottom-left Microhomology",
       title = "Model of mtDNA contacts and deletions") +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank()
  )

cowplot::save_plot(
  plot = pltHeatmap_mhAIC_wt_deletions,
  base_height = 8.316,
  base_width = 11.594,
  file = normalizePath(
    file.path(plots_dir, 'heatmap_microhomology_AIC_wt_deledions.pdf')
  )
)

