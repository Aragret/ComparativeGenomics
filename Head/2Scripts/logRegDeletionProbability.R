mtBreak = read.csv('../../Body/1Raw/MitoBreakDB_12122019.csv')

konstantinData1 = read.table('../../Body/1Raw/compare.square.200.38windows.txt', header=TRUE, sep='\t')
konstantinData2 = read.table('../../Body/1Raw/compare.square.200.38windows_2.txt', header=TRUE, sep='\t')

row.names(konstantinData1) = konstantinData1[, 1]
konstantinData1 = konstantinData1[, -1]

row.names(konstantinData2) = konstantinData2[, 1]
konstantinData2 = konstantinData2[, -1]

a = konstantinData2 / konstantinData1

mtBreak$X5..breakpoint = as.integer(as.character(mtBreak$X5..breakpoint))
mtBreak$X3..breakpoint = as.integer(as.character(mtBreak$X3..breakpoint))

a[is.na(a)] = 0

one_line = c()
for(j in 1:ncol(a)){
  for(i in 1:nrow(a)){
    # i = 3
    # j = 1
    xbegin = as.integer(row.names(a)[i])
    ybegin = as.integer(sub('X', '', colnames(a)[j]))
    xend = as.integer(row.names(a)[i + 1])
    yend = as.integer(sub('X', '', colnames(a)[j + 1]))
    value = a[i,j]
    if(xbegin >= ybegin){
      if(any(mtBreak$X3..breakpoint > xbegin & mtBreak$X3..breakpoint < xend, na.rm = TRUE) &
         any(mtBreak$X5..breakpoint > ybegin & mtBreak$X5..breakpoint < yend, na.rm = TRUE)){
        one_line = rbind(one_line, c(ybegin, xbegin, value, 1))
      }
      else{one_line = rbind(one_line, c(ybegin, xbegin, value, 0))}
    }
  }
}

data = as.data.frame(one_line)
names(data) = c('X5', 'X3', 'Stability', 'Deletions')

for(i in 1:nrow(data)){
  # i = 1
  begin = data[i,]$X5
  end = data[i,]$X3
  if(begin >= 6000 & begin <= 8021 & end >= 14000 & end <= 16021){
    data[i, 'Bublik'] = 1
  }
  else{data[i, 'Bublik'] = 0}
}

data$Deletions = as.factor(data$Deletions)
data$Bublik = as.factor(data$Bublik)

summary(glm(data$Deletions ~ data$Stability + data$Bublik, family = binomial()))

# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)      1.1856     0.1482   7.999 1.25e-15 ***
#   data$Stability  -0.8346     1.0565  -0.790 0.429535    
# data$Bublik1     1.2312     0.3227   3.816 0.000136 *** 


###############################################################################

one_line = c()
for(j in 1:ncol(a)){
  for(i in 1:nrow(a)){
    # i = 3
    # j = 1
    xbegin = as.integer(row.names(a)[i])
    ybegin = as.integer(sub('X', '', colnames(a)[j]))
    xend = as.integer(row.names(a)[i + 1])
    yend = as.integer(sub('X', '', colnames(a)[j + 1]))
    value = a[i,j]
    if(xbegin >= ybegin){
      if(any(mtBreak$X3..breakpoint >= xbegin & mtBreak$X3..breakpoint <= xend &
         mtBreak$X5..breakpoint >= ybegin & mtBreak$X5..breakpoint <= yend, na.rm = TRUE)){
        one_line = rbind(one_line, c(ybegin, xbegin, value, 1))
      }
      else{one_line = rbind(one_line, c(ybegin, xbegin, value, 0))}
    }
  }
}

data = as.data.frame(one_line)
names(data) = c('X5', 'X3', 'Stability', 'Deletions')

summary(data$Deletions)

for(i in 1:nrow(data)){
  # i = 1
  begin = data[i,]$X5
  end = data[i,]$X3
  if(begin >= 6000 & begin <= 8021 & end >= 14000 & end <= 16021){
    data[i, 'Bublik'] = 1
  }
  else{data[i, 'Bublik'] = 0}
}

data$Deletions = as.factor(data$Deletions)
summary(data$Deletions)
# 0    1 
# 1077  354

data$Bublik = as.factor(data$Bublik)

summary(glm(data$Deletions ~ data$Stability + data$Bublik, family = binomial()))

# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)     -1.5583     0.1563  -9.970   <2e-16 ***
#   data$Stability   1.7240     1.0935   1.577    0.115    
# data$Bublik1     2.2633     0.2136  10.594   <2e-16 ***


#########################################################################
# number of deletions

one_line = c()
for(j in 1:ncol(a)){
  for(i in 1:nrow(a)){
    # i = 3
    # j = 1
    xbegin = as.integer(row.names(a)[i])
    ybegin = as.integer(sub('X', '', colnames(a)[j]))
    xend = as.integer(row.names(a)[i + 1])
    yend = as.integer(sub('X', '', colnames(a)[j + 1]))
    value = a[i,j]
    if(xbegin >= ybegin){
      deletionsNumber = sum(mtBreak$X3..breakpoint >= xbegin & mtBreak$X3..breakpoint <= xend &
                              mtBreak$X5..breakpoint >= ybegin & mtBreak$X5..breakpoint <= yend, na.rm = TRUE)
      if(deletionsNumber > 0){
        one_line = rbind(one_line, c(ybegin, xbegin, value, deletionsNumber))
      }
      else{one_line = rbind(one_line, c(ybegin, xbegin, value, 0))}
    }
  }
}

data = as.data.frame(one_line)
names(data) = c('X5', 'X3', 'Stability', 'Deletions')

summary(data$Deletions)

for(i in 1:nrow(data)){
  # i = 1
  begin = data[i,]$X5
  end = data[i,]$X3
  if(begin >= 6000 & begin <= 8021 & end >= 14000 & end <= 16021){
    data[i, 'Bublik'] = 1
  }
  else{data[i, 'Bublik'] = 0}
}

# data$Deletions = as.factor(data$Deletions)
data$Bublik = as.factor(data$Bublik)


summary(lm(data$Deletions ~ data$Stability + data$Bublik))

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      0.4752     0.0869   5.469 5.34e-08 ***
#   data$Stability  -0.2832     0.6266  -0.452    0.651    
# data$Bublik1     1.3568     0.1278  10.616  < 2e-16 ***
