install.packages("SearchTrees")
install.packages("RANN")
install.packages("LS2Wstat")
install.packages("deamer")
# library(SearchTrees)
# library(RANN)
# library(LS2Wstat)
library(deamer) #rlaplace

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


### 

# tree <- list(list(1, 2), list(3, list(4, 5)))
# 
# 
# ##### 1D range query algorithm #####
# 
# findSplitNode <- function(tree, start, end)
# 
# runif(20, 1,100)
# dataset <- sample(1:10,10,replace=F)
# sort(dataset)
# 
# 
# RangeQuery_1D <- function(dataset, range)
# {
#   
# }


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

for(i in m:1){
  for(j in 1:m){
    grids[i,j] <- nrow(testdata[which(
                                      (minx + (j-1)*gridx) < testdata$V1 & 
                                      testdata$V1 < (minx + j*gridx) & 
                                      (miny + (m-i)*gridy) < testdata$V2 & 
                                      testdata$V2 < (miny + (8-i)*gridy)
                                      ),])
  }
}

# for(i in 1:m){
#   for(j in 1:m){
#     grids[i,j] <- nrow(testdata[which((min(testdata$V3) + (j-1)*v/m) < testdata$V3 & testdata$V3 < (min(testdata$V3) + j*v/m) & 
#                                         (min(testdata$V4) + (i-1)*h/m) < testdata$V4 & testdata$V4 < (min(testdata$V4) + i*h/m)),])
#   }
# }


# add noise
noisedgrids <- matrix(0,m,m)
for(i in 1:m){
  for(j in 1:m){
    noisedgrids[i,j] <- grids[i,j] + rlaplace(n = 1, mu = 0, b = 1/e)
  }
}

# 
# grid1 <- testdata[which(testdata$V3 < (min(testdata$V3) + v/m) & 
#                         testdata$V4 < (min(testdata$V4) + h/m)),]
# 
# grid2 <- testdata[which((min(testdata$V3) + 5*v/m) < testdata$V3 & testdata$V3 < (min(testdata$V3) + 6*v/m) & 
#                           (min(testdata$V4) + 0*h/m) < testdata$V4 & testdata$V4 < (min(testdata$V4) + 1*h/m)),]
# 
# plot(grid2, xlim=c(min(testdata[,1]),max(testdata[,1])), 
#             ylim=c(min(testdata[,2]),max(testdata[,2])))

# draw grids
for(i in 0:m){
  abline(v = minx + i*gridx)
  abline(h = miny + i*gridy)
}

# test
xmin = 0
xmax = 60
ymin= -50
ymax = 100

abline(v = xmin, col="red")
abline(v = xmax, col="red")
abline(h = ymin, col="red")
abline(h = ymax, col="red")

# range query original
orq <- function(xmin,xmax,ymin,ymax,dataset){
  results <- dataset[which(dataset$V1 > xmin & dataset$V1 < xmax & 
                             dataset$V2 > ymin & dataset$V2 < ymax),]
  numberOfPoints <- nrow(results)
  return(numberOfPoints)
}

# private range query
prq <- function(xmin,xmax,ymin,ymax,privatedataset){
  numberOfPoints = 0
  i = 1
  while(xmin > minx + (i-1)*gridx){
    i = i + 1
  }
  gridminx = i
  
  i = 1
  while(ymin > miny + (i-1)*gridy){
    i = i + 1
  }
  gridminy = i
  
  j = m
  while(xmax < maxx - (m-j)*gridx){
    j = j - 1
  }
  gridmaxx = j
  
  j = m
  while(ymax < maxy - (m-j)*gridy){
    j = j - 1
  }
  gridmaxy = j
}

points = 0
for(x in gridminx:gridmaxx){
  for(y in gridminy:gridmaxy)
    points = points + grids[x,y] 
}

# upbound
for(x in gridminx:gridmaxx){
  y = gridminy
  numberOfPoints = grids[x,y]*((ymax-gridmaxy*gridy)/gridy)
}


