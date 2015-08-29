##### Packages ####################

install.packages("SearchTrees")
install.packages("RANN")
install.packages("LS2Wstat")
install.packages("deamer")
install.packages("OOmisc")
install.packages("data.tree")
# library(SearchTrees)
# library(RANN)
# library(LS2Wstat)
library(deamer) #rlaplace
library(OOmisc) #rlaplace
# library(data.tree)

##### Data Sampling, Input and Plot ####################

# random sample subset the data set
dataset <- read.table("loc-gowalla_totalCheckins.txt")
s <- sample(nrow(dataset), 100000)
n <- dataset[s,]
n <- n[,3:4]
plot(n)

testdata <- n

write.table(testdata, "testdata100000.txt", sep="\t", row.names = FALSE, col.names = FALSE)

##### Original Range Query ####################

orq <- function(rangeminx,rangemaxx,rangeminy,rangemaxy,dataset){
  
  if(rangeminx < minx){
    rangeminx = minx
  }
  
  if(rangemaxx > maxx){
    rangemaxx = maxx
  }
  
  if(rangeminy < miny){
    rangeminy = miny
  }
  
  if(rangemaxy > maxy){
    rangemaxy = maxy
  }
  
  results <- dataset[which(dataset$V1 >= rangeminx & dataset$V1 <= rangemaxx & 
                             dataset$V2 >= rangeminy & dataset$V2 <= rangemaxy),]
  
  numberOfPoints <- nrow(results)
  
  return(numberOfPoints)
}

#################### UG ####################

testdata <- read.table("testdata100000.txt")

e = 0.9 # epsilon
a = 0.1 # very small portion of the total epsilon
N = nrow(testdata) + rlaplace(n = 1, mu = 0, b = 1/(a*e))  # estimate number of data points
c = 10 # constant(can be changed)

m <- sqrt(N*e/c) # number of grids in both coordinates
m <- round(m,0) # 0 decimal places

m

minx <- min(testdata$V1) #minimum coordinate x
maxx <- max(testdata$V1) #maximum coordinate x
miny <- min(testdata$V2) #minimum coordinate y
maxy <- max(testdata$V2) #maximum coordinate y

plot(testdata, xlim=c(minx,maxx), ylim=c(miny,maxy))

v1 <- maxx - minx
v2 <- maxy - miny

gridx = v1/m #length of each grid
gridy = v2/m #height of each grid

# create grids
grids <- matrix(0,m,m)

for(i in 1:m){
  for(j in 1:m){
    if(i == 1 || j == 1){
      grids[i,j] <- nrow(testdata[which(
        testdata$V1 >= (minx + (i-1)*gridx) & 
          testdata$V1 <= (minx + i*gridx) & 
          testdata$V2 >= (miny + (j-1)*gridy) & 
          testdata$V2 <= (miny + j*gridy)
      ),])
    }else{grids[i,j] <- nrow(testdata[which(
      testdata$V1 > (minx + (i-1)*gridx) & 
        testdata$V1 <= (minx + i*gridx) & 
        testdata$V2 > (miny + (j-1)*gridy) & 
        testdata$V2 <= (miny + j*gridy)
    ),])
    }
  }
  print(i)
}

# add noise only has (1-a)*e privacy baget
noisedgrids <- matrix(0,m,m)
for(i in 1:m){
  for(j in 1:m){
    noisedgrids[i,j] <- grids[i,j] + rlaplace(n = 1, mu = 0, b = 1/((1-a)*e))
  }
}

# draw grids
for(i in 0:m){
  abline(v = minx + i*gridx)
  abline(h = miny + i*gridy)
}

##### Private Range Query For UG ####################

prq <- function(rangeminx,rangemaxx,rangeminy,rangemaxy,privatedataset){
  
  if(rangeminx < minx){
    rangeminx = minx
  }
  
  if(rangemaxx > maxx){
    rangemaxx = maxx
  }
  
  if(rangeminy < miny){
    rangeminy = miny
  }
  
  if(rangemaxy > maxy){
    rangemaxy = maxy
  }
  
  gridminx = 1
  gridminy = 1
  gridmaxx = m
  gridmaxy = m
  
  while(rangeminx > minx + (gridminx-1)*gridx){
    gridminx = gridminx + 1
  }
  
  while(rangeminy > miny + (gridminy-1)*gridy){
    gridminy = gridminy + 1
  }
  
  while(rangemaxx < maxx - (m-gridmaxx)*gridx){
    gridmaxx = gridmaxx - 1
  }
  
  while(rangemaxy < maxy - (m-gridmaxy)*gridy){
    gridmaxy = gridmaxy - 1
  }
  
  if(gridminx > gridmaxx){
    gridminx = gridmaxx
  }
  
  if(gridminy > gridmaxy){
    gridminy = gridmaxy
  }
  
  # topbound
  toppoints = 0
  if(rangemaxy < maxy){
    for(x in gridminx:gridmaxx){
      y = gridmaxy + 1
      toppoints = toppoints + privatedataset[x,y]*((rangemaxy-(miny+gridmaxy*gridy))/gridy)
    }
  }
  
  # downbound
  if(rangeminy > miny){
    downpoints = 0
    for(x in gridminx:gridmaxx){
      y = gridminy - 1
      downpoints = downpoints + privatedataset[x,y]*(((miny+gridminy*gridy)-rangeminy)/gridy)
    }
  }
  
  # left
  leftpoints = 0
  if(rangeminx > minx){
    for(y in gridminy:gridmaxy){
      x = gridminx - 1
      leftpoints = leftpoints + privatedataset[x,y]*(((minx+gridminx*gridx)-rangeminx)/gridx)
    }
  }
  
  # right
  rightpoints = 0
  if(rangemaxx < maxx){
    for(y in gridminy:gridmaxy){
      x = gridmaxx + 1
      rightpoints = rightpoints + privatedataset[x,y]*((rangemaxx - (minx+gridmaxx*gridx))/gridx)
    }
  }
  
  # core
  corepoints = 0
  for(i in gridminx:gridmaxx){
    for(j in gridminy:gridmaxy){
      corepoints = corepoints + privatedataset[i,j]
    }
  }
  
  # tlcorner
  tlpoints = 0
  if(rangeminx > minx && rangemaxy < maxy){
    tlpoints = ((((minx+(gridminx-1)*gridx)-rangeminx)*(rangemaxy-(miny+gridmaxy*gridy)))/(gridx*gridy))*privatedataset[gridminx-1,gridmaxy+1]
  }
  # trcorner
  trpoints = 0
  if(rangemaxx < maxx && rangemaxy < maxy){
    trpoints = (((rangemaxx-(minx+(gridmaxx)*gridx))*(rangemaxy-(miny+(gridmaxy)*gridy)))/(gridx*gridy))*privatedataset[gridmaxx+1,gridmaxy+1]
  }
  # blcorner
  blpoints = 0
  if(rangeminx > minx && rangeminy > miny){
    blpoints = ((((minx+(gridminx-1)*gridx)-rangeminx)*((miny+(gridminy-1)*gridy)-rangeminy))/(gridx*gridy))*privatedataset[gridminx-1,gridminy-1]
  }
  # brcorner
  brpoints = 0
  if(rangemaxx < maxx && rangeminy > miny){
    brpoints = (((rangemaxx-(minx+(gridmaxx)*gridx))*((miny+(gridminy-1)*gridy)-rangeminy))/(gridx*gridy))*privatedataset[gridmaxx+1,gridminy-1]
  }
  
  numberOfPoint = tlpoints+trpoints+blpoints+brpoints+toppoints+downpoints+leftpoints+rightpoints+corepoints
  
  return(numberOfPoint)
}

##### Random Query with spesific area ####################

mkqry <- function(minarea, maxarea, numberOfQries){
  randomquerues <- c()
  while(is.null(nrow(randomquerues))){
    randx <- sort(runif(2,minx,maxx))
    randy <- sort(runif(2,miny,maxy))
    if((randx[2]-randx[1])*(randy[2]-randy[1]) >= minarea &&
      (randx[2]-randx[1])*(randy[2]-randy[1]) <= maxarea){
      randomquerues <- rbind(randomquerues, c(randx,randy))
    }
  }
  
  while(nrow(randomquerues) < numberOfQries){
    randx <- sort(runif(2,minx,maxx))
    randy <- sort(runif(2,miny,maxy))
    if((randx[2]-randx[1])*(randy[2]-randy[1]) >= minarea &&
      (randx[2]-randx[1])*(randy[2]-randy[1]) <= maxarea){
      randomquerues <- rbind(randomquerues, c(randx,randy))
    }
  }
  
  return(randomquerues)
}

totalarea <- (maxx-minx)*(maxy-miny)
maxarea <- totalarea/2
minarea <- totalarea/4

n = 200 # number of queries wanna generate
randboxes <- mkqry(minarea, maxarea, n)

# plot the boxes
plot(testdata, xlim=c(minx,maxx), ylim=c(miny,maxy))

for(i in 1:nrow(randboxes)){
  rect(randboxes[i,1], randboxes[i,3], randboxes[i,2], randboxes[i,4], border = "red")
}

##### RelativeError ####################
relativeError <- c()
for(i in 1:nrow(randboxes)){
  pa = prq(randboxes[i,1], randboxes[i,2], randboxes[i,3], randboxes[i,4],noisedgrids)
  oa = orq(randboxes[i,1], randboxes[i,2], randboxes[i,3], randboxes[i,4],testdata)
  p = 0.001*N
  relativeError <- c(relativeError, (abs(pa-oa))/max(oa,p))
}
length(relativeError[which(relativeError < 10)])
boxplot(relativeError[which(relativeError < 1)])


pa = prq(randcoorminx,randcoormaxx,randcoorminy,randcoormaxy,noisedgrids)
oa = orq(randcoorminx,randcoormaxx,randcoorminy,randcoormaxy,testdata)
p = 0.001*N
relativeError = (abs(pa-oa))/max(oa,p)
relativeError

##### TEST ####################

rangeminx = 0
rangemaxx = 60
rangeminy= 0
rangemaxy = 100

abline(v = rangeminx, col="red")
abline(v = rangemaxx, col="red")
abline(h = rangeminy, col="red")
abline(h = rangemaxy, col="red")

#################### AG #################### 

testdata <- read.table("testdata5000.txt")
sa = 0.01
N = nrow(testdata) + rlaplace(n = 1, mu = 0, b = 1/(sa*e))  # estimate number of data points
a = 0.5 # [0.2,0.6] from the paper
e = 0.1 # epcilon

##### AG level 1 ####################

plot(testdata)
m1 <- max(10,sqrt(N*e/c)/4)


# add noise
noisedgrids <- matrix(0,m1,m1)
for(i in 1:m1){
  for(j in 1:m1){
    noisedgrids[i,j] <- grids[i,j] + rlaplace(n = 1, mu = 0, b = 1/(a*(1-sa)*e))
  }
}

##### AG level 2 ####################

m2 <- matrix(0,m1,m1)
c2 = c/2

for(i in 1:m1){
  for(i in 1:m1){
    m2[i,j] <- sqrt(noisedgrids[i,j]*(1-a)*(1-sa)*e/c2)
  }
}

#################### PSD Based On Full Quadtree #################### 

testdata <- read.table("testdata5000.txt")

minx <- min(testdata$V1) #minimum coordinate x
maxx <- max(testdata$V1) #maximum coordinate x
miny <- min(testdata$V2) #minimum coordinate y
maxy <- max(testdata$V2) #maximum coordinate y

plot(testdata, xlim=c(minx,maxx), ylim=c(miny,maxy))

v1 <- maxx - minx
v2 <- maxy - miny

h = 3 #level of the full quadtree

abline(v=minx, h=miny)
abline(v=maxx, h=maxy)
abline(v=minx+v1/2, h=miny+v2/2,col="red")

##### 4 Part Of The Quadtree #################### 

nw <- function(dataset){
  nwdata <- dataset[which(
    dataset$V1 >= min(dataset$V1) &
      dataset$V1 <= min(dataset$V1) + (max(dataset$V1)-min(dataset$V1))/2 &
      dataset$V2 >= min(dataset$V2) + (max(dataset$V2)-min(dataset$V2))/2 &
      dataset$V2 <= max(dataset$V2)
  ),]
  return(nwdata)
}

ne <- function(dataset){
  nedata <- dataset[which(
    dataset$V1 > min(dataset$V1) + (max(dataset$V1)-min(dataset$V1))/2 &
      dataset$V1 <= max(dataset$V1) &
      dataset$V2 > min(dataset$V2) + (max(dataset$V2)-min(dataset$V2))/2 &
      dataset$V2 <= max(dataset$V2)
  ),]
  return(nedata)
}

sw <- function(dataset){
  swdata <- dataset[which(
    dataset$V1 >= min(dataset$V1) &
      dataset$V1 < min(dataset$V1) + (max(dataset$V1)-min(dataset$V1))/2 &
      dataset$V2 >= min(dataset$V2) &
      dataset$V2 < min(dataset$V2) + (max(dataset$V2)-min(dataset$V2))/2
  ),]
  return(swdata)
}

se <- function(dataset){
  sedata <- dataset[which(
    dataset$V1 >= min(dataset$V1) + (max(dataset$V1)-min(dataset$V1))/2 &
      dataset$V1 <= max(dataset$V1) &
      dataset$V2 >= min(dataset$V2) &
      dataset$V2 <= min(dataset$V2) + (max(dataset$V2)-min(dataset$V2))/2
  ),]
  return(sedata)
}

##### Data Decomposition Based On Full Quadtree #################### 

quad <- function(dataset, h){
  quaddata <- matrix(list(), nrow = h, ncol=4^(h-1))
  quaddata[[1,1]] <- dataset
  for(i in 2:h){
    for(j in seq(1, 4^(i-1), by=4)){
      quaddata[[i,j]] <- nw(quaddata[[i-1,(j+3)/4]])
      quaddata[[i,j+1]] <- ne(quaddata[[i-1,(j+3)/4]])
      quaddata[[i,j+2]] <- sw(quaddata[[i-1,(j+3)/4]])
      quaddata[[i,j+3]] <- se(quaddata[[i-1,(j+3)/4]])
    }
    print(j)
  }
  return(quaddata)
}

##### Generate Noisetree (only has count) #################### 

noisetree <- function(dataset, h, te){
  # te - total epsilon
  # h - height of the tree 
  # choose epsilon
  e <- c()
  for(i in 1:h){
    # e[i] <- (2^((h-i)/3))*te*((2^(1/3)-1)/(2^((h+1)/3)-1))
    e[i] <- te/h
  }
  
  # creat noise tree
  noisetree <- matrix(0, nrow = h, ncol=4^(h-1))
  tree <- quad(dataset, h)
  
  noisetree[1,1] <- nrow(tree[[1,1]]) + rlaplace(1, mu=0, b=1/e[1])
  
  for(i in 2:h){
    for(j in 1:4^(i-1)){
      noisetree[i,j] <- nrow(tree[[i,j]]) + rlaplace(1, mu=0, b=1/e[i])
    }
  }
  
  return(noisetree)
}

##### Private Range Query For PSD ####################

prqpsd <- function(rangeminx,rangemaxx,rangeminy,rangemaxy,dataset){
  
  if(rangeminx < min(dataset$V1)){
    rangeminx = min(dataset$V1)
  }
  
  if(rangemaxx > max(dataset$V1)){
    rangemaxx = max(dataset$V1)
  }
  
  if(rangeminy < min(dataset$V2)){
    rangeminy = min(dataset$V2)
  }
  
  if(rangemaxy > max(dataset$V2)){
    rangemaxy = max(dataset$V2)
  }

  ft <- fullquadrange(dataset, h)
  noisecounts = 0
  minpt <- c()
  maxpt <- c()
  for(i in h){
    for(j in 1:4^(i-1)){
      if(rangeminx >= ft[[i,j]][1,1] && rangeminy >= ft[[i,j]][1,2]&&
         rangeminx <= ft[[i,j]][2,1] && rangeminy <= ft[[i,j]][2,2]){
        minpt <- rbind(minpt, c(i,j))
      }
      
      if(rangemaxx >= ft[[i,j]][1,1] && rangemaxy >= ft[[i,j]][1,2]&&
         rangemaxx <= ft[[i,j]][2,1] && rangemaxy <= ft[[i,j]][2,2]){
        maxpt <- rbind(maxpt, c(i,j))
      }
    }
  }
  from <- ft[minpt][[1]][1,]
  to <- ft[maxpt][[1]][2,]
  
########## all grids that containd and intersects with range query
  allgridsinvloved <- c()
  for(j in 1:4^(h-1)){
    if(from[1] <= ft[[h,j]][1,1] && from[2] <= ft[[h,j]][1,2] &&
       to[1] >= ft[[h,j]][2,1] && to[2] >= ft[[h,j]][2,2]){
      allgridsinvloved <- rbind(allgridsinvloved, c(h,j))
    }
  }

########## fully contained grids with minimun noised count (optimun path + any other ones)

  x <- c()
  for(i in 1:h){
    for(j in 1:4^(i-1)){
      if(ft[[i,j]][1,1] >= rangeminx && ft[[i,j]][1,2] >= rangeminy &&
         ft[[i,j]][2,1] <= rangemaxx && ft[[i,j]][2,2] <= rangemaxy){
        x <- rbind(x, c(i,j))
      }
    }
  }
  
########## edges grid
  
  fullyunitgrids <- x[which(x[,1] == h),]
  if(!is.matrix(fullyunitgrids)){
    fullyunitgrids <- t(as.matrix(fullyunitgrids))
  }
  newx <- c()
  
  for(i in 1:nrow(allgridsinvloved)){
    for(j in 1:nrow(fullyunitgrids)){
      if(allgridsinvloved[i,1] == fullyunitgrids[j,1] && allgridsinvloved[i,2] == fullyunitgrids[j,2]){
        newx <- c(newx, i)
      }
    }
  }
  edgegrids <- allgridsinvloved[-newx,]
  
########## corners grid
  
  bl <- c()
  br <- c()
  tl <- c()
  tr <- c()
  cornergrids <- c()
  
  for(j in 1:4^(h-1)){
    if(ft[[h,j]][1,1] < rangeminx && ft[[h,j]][1,2] < rangeminy &&
       ft[[h,j]][2,1] > rangeminx && ft[[h,j]][2,2] > rangeminy){
      bl <- c(h,j)
      bl <- t(as.matrix(bl))
      blfrac = (ft[[h,j]][2,1]-rangeminx)*(ft[[h,j]][2,2]-rangeminy)/(gridv1*gridv2)
      blnoisecount = blfrac*nt[h,j]
    }
    if(ft[[h,j]][1,1] < rangemaxx && ft[[h,j]][1,2] < rangeminy &&
       ft[[h,j]][2,1] > rangemaxx && ft[[h,j]][2,2] > rangeminy){
      br <- c(h,j)
      br <- t(as.matrix(br))
      brfrac = (rangemaxx-ft[[h,j]][1,1])*(rangeminy-ft[[h,j]][1,2])/(gridv1*gridv2)
      brnoisecount = brfrac*nt[h,j]
    }
    if(ft[[h,j]][1,1] < rangeminx && ft[[h,j]][1,2] < rangemaxy &&
       ft[[h,j]][2,1] > rangeminx && ft[[h,j]][2,2] > rangemaxy){
      tl <- c(h,j)
      tl <- t(as.matrix(tl))
      tlfrac = (ft[[h,j]][2,1]-rangeminx)*(rangemaxy-ft[[h,j]][1,2])/(gridv1*gridv2)
      tlnoisecount = tlfrac*nt[h,j]
    }
    if(ft[[h,j]][1,1] < rangemaxx && ft[[h,j]][1,2] < rangemaxy &&
       ft[[h,j]][2,1] > rangemaxx && ft[[h,j]][2,2] > rangemaxy){
      tr <- c(h,j)
      tr <- t(as.matrix(tr))
      trfrac = (rangemaxx-ft[[h,j]][1,1])*(rangemaxy-ft[[h,j]][1,2])/(gridv1*gridv2)
      trnoisecount = trfrac*nt[h,j]
    }
  }
  
  cornergrids <- rbind(tl,tr,bl,br)
  noisedcorners = trnoisecount+tlnoisecount+brnoisecount+blnoisecount
  
########## top grids
  
  topgrids <- c()
  topnoisecount = 0
  for(j in 1:4^(h-1)){
    if(ft[[h,j]][1,2] < rangemaxy && ft[[h,j]][2,2] > rangemaxy &&
       ft[[h,j]][1,1] > rangeminx && ft[[h,j]][2,1] < rangemaxx){
      topgrids <- rbind(topgrids, c(h,j))
      topnoisecount = topnoisecount+((rangemaxy-ft[[h,j]][1,2])/gridv2)*nt[h,j]
    }
  }
  
########## bottum grids
  
  botgrids <- c()
  botnoisecount = 0
  for(j in 1:4^(h-1)){
    if(ft[[h,j]][1,2] < rangeminy && ft[[h,j]][2,2] > rangeminy &&
       ft[[h,j]][1,1] > rangeminx && ft[[h,j]][2,1] < rangemaxx){
      botgrids <- rbind(botgrids, c(h,j)) 
      botnoisecount = botnoisecount+((ft[[h,j]][2,2]-rangeminy)/gridv2)*nt[h,j]
    }
  }

########## left grids
  
  leftgrids <- c()
  leftnoisecount = 0
  for(j in 1:4^(h-1)){
    if(ft[[h,j]][1,2] > rangeminy && ft[[h,j]][2,2] < rangemaxy &&
       ft[[h,j]][1,1] < rangeminx && ft[[h,j]][2,1] > rangeminx){
      leftgrids <- rbind(leftgrids, c(h,j))
      leftnoisecount = leftnoisecount+((ft[[h,j]][2,1]-rangeminx)/gridv1)*nt[h,j]
    }
  }
  
########## right grids
  
  rightgrids <- c()
  rightnoisecount = 0
  for(j in 1:4^(h-1)){
    if(ft[[h,j]][1,2] > rangeminy && ft[[h,j]][2,2] < rangemaxy &&
       ft[[h,j]][1,1] < rangemaxx && ft[[h,j]][2,1] > rangemaxx){
      rightgrids <- rbind(rightgrids, c(h,j))
      rightnoisecount = rightnoisecount+((rangemaxx-ft[[h,j]][1,1])/gridv1)*nt[h,j]
    }
  }
  
  foursidenoisecount = rightnoisecount+leftnoisecount+topnoisecount+botnoisecount
  
########## fully contained grids with minimun noised count (ONLY optimun path)
  
  if(nrow(x) >= 4){
    tmp <- x[which(x[,1] != h),]
    if(!is.matrix(tmp)){
      tmp <- t(as.matrix(tmp))
    }
    parent <- tmp[,2]

    n <- length(parent)
    
    for(l in h:1){
      cl <- x[which(x[,1] == l),]
      if(!is.matrix(cl)){
        cl <- t(as.matrix(cl))
      }
      rmchild <- c()
      
      for(nodes in seq(1, 4^(h-1), by = 4)){
        if((c(nodes, nodes+1, nodes+2,nodes+3) %in% cl[,2])[1] &&
           (c(nodes, nodes+1, nodes+2,nodes+3) %in% cl[,2])[2] &&
           (c(nodes, nodes+1, nodes+2,nodes+3) %in% cl[,2])[3] &&
           (c(nodes, nodes+1, nodes+2,nodes+3) %in% cl[,2])[4]){
          rmchild <- rbind(rmchild, cl[which(cl[,2] == nodes:(nodes+3)),])
        }
      }
      
      if(!is.null(nrow(rmchild))){
        newx <- c()
        for(i in 1:nrow(x)){
          for(j in 1:nrow(rmchild)){
            if(x[i,1] == rmchild[j,1] && x[i,2] == rmchild[j,2]){
              newx <- c(newx, i)
            }
          }
        }
        x <- x[-newx,]
      }
    }
  }
  
##########  sum up the noise count for fully contained grids 
 
  for(i in 1:nrow(x)){
    noisecounts = noisecounts + nt[x[i,1],x[i,2]]
  }
  
  totalnoisecounts = noisecounts+noisedcorners+foursidenoisecount
  
  return(totalnoisecounts)
}

##### get range of all grids at each level ####################

minx <- min(dataset$V1)
maxx <- max(dataset$V1)
miny <- min(dataset$V2)
maxy <- max(dataset$V2)

quadrange <- function(minx, miny, maxx, maxy){
  
  datarange <- rbind(c(minx, miny), c(maxx,maxy))
  
  v1 <- datarange[2,1]-datarange[1,1]
  v2 <- datarange[2,2]-datarange[1,2]
  
  nwrange <- rbind(c(datarange[1,1], datarange[1,2]+v2/2), 
                   c(datarange[1,1]+v1/2,datarange[2,2]))
  
  nerange <- rbind(c(datarange[1,1]+v1/2, datarange[1,2]+v2/2), 
                   c(datarange[2,1], datarange[2,2]))
  
  swrange <- rbind(c(datarange[1,1], datarange[1,2]), 
                   c(datarange[1,1]+v1/2, datarange[1,2]+v2/2))
  
  serange <- rbind(c(datarange[1,1]+v1/2, datarange[1,2]), 
                   c(datarange[2,1], datarange[1,2]+v2/2))
  
  range <- list(nwrange, nerange, swrange, serange)
  
  return(range)
}

fullquadrange <- function(dataset, h){
  
  minx <- min(dataset$V1)
  maxx <- max(dataset$V1)
  miny <- min(dataset$V2)
  maxy <- max(dataset$V2)
  
  range <- matrix(list(), nrow = h, ncol = 4^(h-1))
  range[[1,1]] <- rbind(c(minx, miny), c(maxx, maxy))
  
  for(i in 2:h){
    for(j in seq(1, 4^(i-1), by = 4)){
      range[[i,j]] <- quadrange(range[[i-1,(j+3)/4]][1,1], range[[i-1,(j+3)/4]][1,2], range[[i-1,(j+3)/4]][2,1], range[[i-1,(j+3)/4]][2,2])[[1]]
      range[[i,j+1]] <- quadrange(range[[i-1,(j+3)/4]][1,1], range[[i-1,(j+3)/4]][1,2], range[[i-1,(j+3)/4]][2,1], range[[i-1,(j+3)/4]][2,2])[[2]]
      range[[i,j+2]] <- quadrange(range[[i-1,(j+3)/4]][1,1], range[[i-1,(j+3)/4]][1,2], range[[i-1,(j+3)/4]][2,1], range[[i-1,(j+3)/4]][2,2])[[3]]
      range[[i,j+3]] <- quadrange(range[[i-1,(j+3)/4]][1,1], range[[i-1,(j+3)/4]][1,2], range[[i-1,(j+3)/4]][2,1], range[[i-1,(j+3)/4]][2,2])[[4]]
    }
  }
  return(range)
}

##### TEST ####################

rangeminx = 0
rangeminy = 0
rangemaxx = 60
rangemaxy = 100

abline(v = rangeminx, col="red")
abline(v = rangemaxx, col="red")
abline(h = rangeminy, col="red")
abline(h = rangemaxy, col="red")

h = 6
nt <- noisetree(testdata, h, 0.5)
m <- sqrt(4^(h-1))
gridv1 <- (max(testdata$V1) - min(testdata$V1))/m
gridv2 <- (max(testdata$V2) - min(testdata$V2))/m

prqpsd(rangeminx,rangemaxx,rangeminy,rangemaxy,testdata) #private
orq(rangeminx,rangemaxx,rangeminy,rangemaxy,testdata) #original