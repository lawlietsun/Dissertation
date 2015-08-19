install.packages("SearchTrees")
install.packages("RANN")
install.packages("LS2Wstat")
install.packages("deamer")
install.packages("OOmisc")
# library(SearchTrees)
# library(RANN)
# library(LS2Wstat)
library(deamer) #rlaplace
library(OOmisc) #rlaplace

##############################
# x1 <- runif(10, 0, 2*pi)
# x2 <- runif(10, 0,3)
# x <- c(1,2,4,3,5,6,m,8,m,9)
# y <- c(9,3,1,m,4,8,2,8,9,6)
# d <- data.frame(x1, x2)
# d <- data.frame(x,y)
# plot(d)
# abline(v = median(d$x))
# d1 <- d[which(d$x > median(d$x)),]
# d2 <- d[which(d$x <= median(d$x)),]
# abline(h = median(d1$y))
# abline(h = median(d2$y))
# 
# d = cbind(x,y)
# tree = createTree(d)
# inrect = rectLookup(tree, xlim = c(0,5), ylim=c(0, 5))
# 
# 
# nearest <- nn2(DATA,DATA)

# BinaryTree
BinSearch <- function(A, value, low, high) {
  if ( high < low ) {
    return(NULL)
  } else {
    mid <- floor((low + high) / 2)
    if ( A[mid] > value )
      BinSearch(A, value, low, mid-1)
    else if ( A[mid] < value )
      BinSearch(A, value, mid+1, high)
    else
      mid
  }
}
##############################

# random sample subset the data set
dataset <- read.table("loc-gowalla_totalCheckins.txt")
s <- sample(nrow(dataset), 100000)
n <- dataset[s,]
n <- n[,3:4]
plot(n)

testdata <- n

write.table(testdata, "testdata100000.txt", sep="\t", row.names = FALSE, col.names = FALSE)

# start ### UG #######################################################################
testdata <- read.table("testdata5000.txt")

# n=1000
# x <- rchisq(n,3)
# b=0.4 # 1/e 
# e <- rlaplace(n, 0, b)
# y <- x + e  #noisy observations with laplace noise
# 
# noise = rlaplace(1, 0, 1/0.5)


######################
e = 0.1 # epsilon
a = 0.01 # very small portion of the total epsilon
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
    grids[i,j] <- nrow(testdata[which(
      testdata$V1 > (minx + (i-1)*gridx) & 
        testdata$V1 < (minx + i*gridx) & 
        testdata$V2 > (miny + (j-1)*gridy) & 
        testdata$V2 < (miny + j*gridy)
    ),])
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

# original range query ##################
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
  
  results <- dataset[which(dataset$V1 > rangeminx & dataset$V1 < rangemaxx & 
                             dataset$V2 > rangeminy & dataset$V2 < rangemaxy),]
  numberOfPoints <- nrow(results)
  return(numberOfPoints)
}

# private range query ################
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
  if(rangemaxy < maxy){
    toppoints = 0
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
  if(rangeminx > minx){
    leftpoints = 0
    for(y in gridminy:gridmaxy){
      x = gridminx - 1
      leftpoints = leftpoints + privatedataset[x,y]*(((minx+gridminx*gridx)-rangeminx)/gridx)
    }
  }
  
  # right
  if(rangemaxx < maxx){
    rightpoints = 0
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

# random query 1 D/128 ~ D/64
x <- runif(1, 0, N/128)
y <- runif(1, 0, N/128)
randcoorminx <- runif(1,minx,maxx)
randcoorminy <- runif(1,miny,maxy)
randcoormaxx <- randcoorminx + x
randcoormaxy <- randcoorminy + y

# relativeError
pa = prq(randcoorminx,randcoormaxx,randcoorminy,randcoormaxy,noisedgrids)
oa = orq(randcoorminx,randcoormaxx,randcoorminy,randcoormaxy,testdata)
p = 0.001*N
relativeError = (abs(pa-oa))/max(oa,p)
relativeError

privatedataset = noisedgrids

randcoorminx = rangeminx
randcoormaxx = rangemaxx
randcoorminy = rangeminy
randcoormaxy = rangemaxy


# test
rangeminx = 0
rangemaxx = 60
rangeminy= -50
rangemaxy = 100

abline(v = rangeminx, col="red")
abline(v = rangemaxx, col="red")
abline(h = rangeminy, col="red")
abline(h = rangemaxy, col="red")


# AG ###############################################

testdata <- read.table("testdata5000.txt")
sa = 0.01
N = nrow(testdata) + rlaplace(n = 1, mu = 0, b = 1/(sa*e))  # estimate number of data points
a = 0.5 # [0.2,0.6] from the paper
e = 0.1 # epcilon

# level 1 ##########

plot(testdata)
m1 <- max(10,sqrt(N*e/c)/4)


# add noise
noisedgrids <- matrix(0,m1,m1)
for(i in 1:m1){
  for(j in 1:m1){
    noisedgrids[i,j] <- grids[i,j] + rlaplace(n = 1, mu = 0, b = 1/(a*(1-sa)*e))
  }
}

# level 2 ##########
m2 <- matrix(0,m1,m1)
c2 = c/2

for(i in 1:m1){
  for(i in 1:m1){
    m2[i,j] <- sqrt(noisedgrids[i,j]*(1-a)*(1-sa)*e/c2)
  }
}
