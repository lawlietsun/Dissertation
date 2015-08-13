install.packages("SearchTrees")
install.packages("RANN")
install.packages("LS2Wstat")
install.packages("deamer")
library(SearchTrees)
library(RANN)
library(LS2Wstat)
library(deamer)

x1 <- runif(10, 0, 2*pi)
x2 <- runif(10, 0,3)
x <- c(1,2,4,3,5,6,7,8,7,9)
y <- c(9,3,1,7,4,8,2,8,9,6)
d <- data.frame(x1, x2)
d <- data.frame(x,y)
plot(d)
abline(v = median(d$x))
d1 <- d[which(d$x > median(d$x)),]
d2 <- d[which(d$x <= median(d$x)),]
abline(h = median(d1$y))
abline(h = median(d2$y))

d = cbind(x,y)
tree = createTree(d)
inrect = rectLookup(tree, xlim = c(0,5), ylim=c(0, 5))


nearest <- nn2(DATA,DATA)

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

tree <- list(list(1, 2), list(3, list(4, 5)))


##### 1D range query algorithm #####

findSplitNode <- function(tree, start, end)

runif(20, 1,100)
dataset <- sample(1:10,10,replace=F)
sort(dataset)


RangeQuery_1D <- function(dataset, range)
{
  
}


# random sample subset the data set
s <- sample(nrow(dataset), 5000)
n <- dataset[s,]
n <- n[,3:4]
plot(n)

testdata <- n

write.table(testdata, "testdata.txt", sep="\t", row.names = FALSE, col.names = FALSE)

n=1000
x <- rchisq(n,3)
b=0.4 # 1/e 
e <- rlaplace(n, 0, b)
y <- x + e  #noisy observations with laplace noise

noise = rlaplace(1, 0, 1/0.5)


##################################################################
e = 0.1 # epsilon
N = nrow(testdata) # number of data points
c = 10 # constant(can be changed)

m <- sqrt(N*e/c) # grid size

plot(testdata, xlim=c(min(testdata$V3),max(testdata$V3)), 
               ylim=c(min(testdata$V4),max(testdata$V4)))

v <- max(testdata$V3) - min(testdata$V3)
h <- max(testdata$V4) - min(testdata$V4)

grids <- matrix(0,7,7)

for(i in 7:1){
  for(j in 1:7){
    grids[i,j] <- nrow(testdata[which((min(testdata$V3) + (j-1)*v/7) < testdata$V3 & testdata$V3 < (min(testdata$V3) + j*v/7) & 
                                   (min(testdata$V4) + (7-i)*h/7) < testdata$V4 & testdata$V4 < (min(testdata$V4) + (8-i)*h/7)),])
  }
}

grid1 <- testdata[which(testdata$V3 < (min(testdata$V3) + v/7) & 
                        testdata$V4 < (min(testdata$V4) + h/7)),]

grid2 <- testdata[which((min(testdata$V3) + 5*v/7) < testdata$V3 & testdata$V3 < (min(testdata$V3) + 6*v/7) & 
                          (min(testdata$V4) + 0*h/7) < testdata$V4 & testdata$V4 < (min(testdata$V4) + 1*h/7)),]

plot(grid2, xlim=c(min(testdata[,1]),max(testdata[,1])), 
            ylim=c(min(testdata[,2]),max(testdata[,2])))

# draw grids
for(i in 0:7){
  abline(h = min(testdata$V4) + i*h/7)
  abline(v = min(testdata$V3) + i*v/7)
}


