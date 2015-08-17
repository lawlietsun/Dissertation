install.packages("SearchTrees")
install.packages("RANN")
install.packages("LS2Wstat")
install.packages("deamer")
# library(SearchTrees)
# library(RANN)
# library(LS2Wstat)
library(deamer) #rlaplace

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
s <- sample(nrow(dataset), 5000)
n <- dataset[s,]
n <- n[,3:4]
plot(n)

testdata <- n

write.table(testdata, "testdata.txt", sep="\t", row.names = FALSE, col.names = FALSE)

# start ##########################################################################
testdata <- read.table("testdata.txt")

# n=1000
# x <- rchisq(n,3)
# b=0.4 # 1/e 
# e <- rlaplace(n, 0, b)
# y <- x + e  #noisy observations with laplace noise
# 
# noise = rlaplace(1, 0, 1/0.5)


######################
e = 0.1 # epsilon
N = nrow(testdata) # number of data points
c = 10 # constant(can be changed)

m <- sqrt(N*e/c) # number of grids in both coordinates
m <- round(m,0) # 0 decimal places

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
}

# add noise
noisedgrids <- matrix(0,m,m)
for(i in 1:m){
  for(j in 1:m){
    noisedgrids[i,j] <- grids[i,j] + rlaplace(n = 1, mu = 0, b = 1/e)
  }
}

# draw grids
for(i in 0:m){
  abline(v = minx + i*gridx)
  abline(h = miny + i*gridy)
}

# test
rangeminx = 0
rangemaxx = 60
rangeminy= -50
rangemaxy = 100

abline(v = rangeminx, col="red")
abline(v = rangemaxx, col="red")
abline(h = rangeminy, col="red")
abline(h = rangemaxy, col="red")

# original range query
orq <- function(rangeminx,rangemaxx,rangeminy,rangemaxy,dataset){
  results <- dataset[which(dataset$V1 > rangeminx & dataset$V1 < rangemaxx & 
                             dataset$V2 > rangeminy & dataset$V2 < rangemaxy),]
  numberOfPoints <- nrow(results)
  return(numberOfPoints)
}

# private range query
prq <- function(rangeminx,rangemaxx,rangeminy,rangemaxy,privatedataset){
  
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

  # topbound
  toppoints = 0
  for(x in gridminx:gridmaxx){
    y = gridmaxy + 1
    toppoints = toppoints + grids[x,y]*((rangemaxy-(miny+gridmaxy*gridy))/gridy)
  }
  
  # downbound
  downpoints = 0
  for(x in gridminx:gridmaxx){
    y = gridminy - 1
    downpoints = downpoints + grids[x,y]*(((miny+gridminy*gridy)-rangeminy)/gridy)
  }
  
  # left
  leftpoints = 0
  for(y in gridminy:gridmaxy){
    x = gridminx - 1
    leftpoints = leftpoints + grids[x,y]*(((minx+gridminx*gridx)-rangeminx)/gridx)
  }
  
  # right
  rightpoints = 0
  for(y in gridminy:gridmaxy){
    x = gridmaxx + 1
    rightpoints = rightpoints + grids[x,y]*((rangemaxx - (minx+gridmaxx*gridx))/gridx)
  }
  
  # core
  corepoints = 0
  for(i in gridminx:gridmaxx){
    for(j in gridminy:gridmaxy){
      corepoints = corepoints + grids[i,j]
    }
  }
  
  # tlcorner
  tlpoints = (((minx+gridminx*gridx)-rangeminx)*(rangemaxy-(miny+gridmaxy*gridy)))/(gridx*gridy)
  # trcorner
  trpoints = ((rangemaxx-(minx+gridmaxx*gridx))*(rangemaxy-(miny+gridmaxy*gridy)))/(gridx*gridy)
  # blcorner
  blpoints = (((minx+gridminx*gridx)-rangeminx)*((miny+gridmaxy*gridy)-rangeminy))/(gridx*gridy)
  # brcorner
  brpoints = ((rangemaxx-(minx+gridmaxx*gridx))*((miny+gridmaxy*gridy)-rangeminy))/(gridx*gridy)
  
  numberOfPoint = tlpoints+trpoints+blpoints+brpoints+toppoints+downpoints+leftpoints+rightpoints+corepoints
  
  return(numberOfPoint)
}

pa = prq(rangeminx,rangemaxx,rangeminy,rangemaxy,privatedataset)
oa = orq(rangeminx,rangemaxx,rangeminy,rangemaxy,testdata)
p = 0.001*N
relativeError = (abs(pa-oa))/max(oa,p)

ee