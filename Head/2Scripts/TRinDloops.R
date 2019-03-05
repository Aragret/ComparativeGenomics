rm(list = ls())

tr = read.table('../../Body/3Results/TandRepInfo.txt', header=TRUE, sep = '\t')
dloops = read.table('../../Body/2Derived/dloops_control_regions.txt', header=TRUE, sep='\t',
                    row.names = NULL)

names(dloops) = names(dloops)[-1]
dloops = dloops[,-7]


######################################################################################
##########################     get start end for each dloop

# a = '[12687:14545](+)'

# coords = c()
# for(i in strsplit(a, ':')){
#   coords = c(coords, gsub('\\D', '', i))
# }

dloops$dlStart = c(1:3619)
dloops$dlEnd = c(1:3619)

for(i in 1:nrow(dloops)){
  # i = 267
  a = as.character(dloops$Feature_location[i])
  coords = c()
  b = strsplit(a, ',')
  
  if(length(b[[1]]) == 1){
    for(j in strsplit(a, ':')){
      coords = c(coords, gsub('\\D', '', j))
    }
    Start = as.numeric(coords[1]); End = as.numeric(coords[2])
  }
  if(length(b[[1]]) == 2){
    first_pare = strsplit(b[[1]][1], ':')
    first = as.numeric(gsub('\\D', '', first_pare[[1]][1]))
    if(first != 0){
      Start = first;
      second_pare = strsplit(b[[1]][2], ':')
      End = as.numeric(gsub('\\D', '', second_pare[[1]][2]))
    }
    if(first == 0){
      End = as.numeric(gsub('\\D', '', first_pare[[1]][2]))
      second_pare = strsplit(b[[1]][2], ':')
      Start = as.numeric(gsub('\\D', '', second_pare[[1]][1]))
    }
  }
  dloops[i, 'dlStart'] = Start; dloops[i, 'dlEnd'] = End
}

######################################################################################
#########################      location of TRs

coords = merge(tr, dloops[, c("Species", "dlStart", "dlEnd")], by='Species', all.x = TRUE)
coords = coords[!is.na(coords$dlStart),]

for(i in 1:nrow(coords)){
  # i = 1
  # 
  if(coords$dlStart[i] < coords$dlEnd[i]){
    if((coords$Start[i] >= coords$dlStart[i]) & (coords$End[i] <= coords$dlEnd[i])){
      coords[i, 'InDloop'] = 1
    }
    if((coords$Start[i] < coords$dlStart[i]) & (coords$End[i] <= coords$dlStart[i]) | 
       (coords$Start[i] >= coords$dlEnd[i]) & (coords$End[i] > coords$dlEnd[i])){
      coords[i, 'InDloop'] = 0
    }
    if((coords$Start[i] < coords$dlStart[i]) & (coords$End[i] > coords$dlStart[i]) | 
       (coords$Start[i] < coords$dlEnd[i]) & (coords$End[i] > coords$dlEnd[i])){
      coords[i, 'InDloop'] = 2 # partially
    }
  }
  
  if(coords$dlStart[i] > coords$dlEnd[i]){
    if((coords$Start[i] < coords$dlEnd[i]) & (coords$End[i] <= coords$dlEnd[i]) |
       (coords$Start[i] >= coords$dlStart[i]) & (coords$End[i] > coords$dlStart[i])){
      coords[i, 'InDloop'] = 1
    }
    if((coords$Start[i] < coords$dlStart[i]) & (coords$End[i] > coords$dlEnd[i])){
      coords[i, 'InDloop'] = 0
    }
    if((coords$Start[i] < coords$dlEnd[i]) & (coords$End[i] > coords$dlEnd[i]) | 
       (coords$Start[i] < coords$dlStart[i]) & (coords$End[i] > coords$dlStart[i])){
      coords[i, 'InDloop'] = 2 # partially
    }
  }
}

table(coords$InDloop)

# a = coords[, c('Start', 'End', 'dlStart', 'dlEnd', 'InDloop')]

write.table(coords[, -c(27, 28)], '../../Body/3Results/TandRepInfo.txt', sep='\t', quote = F,
            row.names = F)
