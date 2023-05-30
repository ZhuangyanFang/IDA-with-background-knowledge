
rm(list=ls())

library(GeneNet)
library(bnlearn)
library(igraph)
library(graph)
library(pcalg)

set.seed(10)
pth <- 'E:\\Research\\bgIDA\\code\\code'
setwd(pth)
source(paste(pth, 'bgIDA.R', sep = "\\", collapse = NULL))

data("arth800")
summary(arth800.expr)

#' # Compute Partial Correlations and Select Relevant Edges
colnames(arth800.expr) <- 1:ncol(arth800.expr)
pcor.dyn = ggm.estimate.pcor(arth800.expr, method = "dynamic")
arth.edges = network.test.edges(pcor.dyn, direct=FALSE)
dim(arth.edges)

# arth.net = extract.network(arth.edges, method.ggm="number", cutoff.ggm=150)

arth.net = extract.network(arth.edges, cutoff.ggm=0.999)
dim(arth.net)

# learning
node.labels = as.character(1:ncol(arth800.expr))
arth.dataframe <- data.frame(matrix(arth800.expr, nrow = nrow(arth800.expr), ncol = ncol(arth800.expr)))
colnames(arth.dataframe) <- colnames(arth800.expr)
fixedGaps = matrix(TRUE, ncol(arth800.expr), ncol(arth800.expr))
colnames(fixedGaps) <- rownames(fixedGaps) <- colnames(arth800.expr)

for (i in 1:nrow(arth.net)){
  x = arth.net[i, 'node1']
  y = arth.net[i, 'node2']
  fixedGaps[x, y] <- fixedGaps[y, x] <- FALSE
}
init.graph <- fixedGaps == FALSE
pos <- which(colSums(init.graph) == 0 & colSums(t(init.graph)) == 0)
init.graph <- init.graph[-pos, -pos]

# setting blacklist
tmp <- which(init.graph == FALSE, arr.ind = T)
black <- matrix(0, nrow(tmp), 2)
select.name <- colnames(init.graph)
black[, 1] <- select.name[tmp[, 1]]
black[, 2] <- select.name[tmp[, 2]]
black <- data.frame(black)

# setting white list, white list are edges from bnlearn
load(paste(pth, 'arth150.rda', sep = "\\", collapse = NULL))
adj.true <- amat(bn)
dag.true <- as(adj.true, 'graphNEL')
cpdag.true <- dag2cpdag(dag.true)
names.true <- colnames(adj.true)
edgeList <- which(adj.true != 0, arr.ind = T)
x <- edgeList[, 1]
y <- edgeList[, 2]
edgeList[, 1] <- names.true[x]
edgeList[, 2] <- names.true[y]

# checking true names and select names
length(intersect(names.true, select.name)) == length(names.true)

# from blacklist remove those edges in the whitelist
pos <- c()
for (i in 1:nrow(edgeList)){
  if (sum(black[, 1] == edgeList[i, 1] & black[, 2] == edgeList[i, 2])>0){
    pos <- c(pos, i)
  }
}
if (length(pos) > 0){
  black <- black[-pos, ]
}
  
white <- data.frame(edgeList)

subdata <- arth.dataframe[, select.name]
res.bn <- tabu(subdata, blacklist = black, whitelist = white, debug = FALSE)
res.adj <- amat(res.bn)
pos <- which(colSums(res.adj) == 0 & colSums(t(res.adj)) == 0)
subsubdata <- subdata[, -pos]
res.dag.adj <- res.adj[-pos, -pos]
res.dag <- as(res.dag.adj, 'graphNEL')
res.cpdag <- dag2cpdag(res.dag)
res.cpdag.adj <- as(res.cpdag, 'matrix')
cat('number of undirected edge:', nrow(which(res.cpdag.adj == 1 & t(res.cpdag.adj) == 1, arr.ind = T))/2, '\n')
cat('number of directed edge:', nrow(which(res.dag.adj == 1, arr.ind = T) )
    - nrow(which(res.cpdag.adj == 1 & t(res.cpdag.adj) == 1, arr.ind = T))/2, '\n')
# checking white list
for (i in 1:nrow(white)){
  if (!(res.dag.adj[as.character(white[i, 1]), as.character(white[i, 2])] == 1 
        & res.dag.adj[as.character(white[i, 2]), as.character(white[i, 1])] == 0)){
    print(i)
    break
  }
} 
  
globalAttrs = list()
globalAttrs$edge = list(color = "black", lty = "solid", lwd = 1, arrowsize=1)
globalAttrs$node = list(fillcolor = "lightblue", shape = "ellipse", fixedsize = FALSE)
plot(res.cpdag, attrs = globalAttrs)


# getting direct causal information
# getting non-ancestral information.

xx <- c()
yy <- c()
anxx <- c()
anyy <- c()
for (i in 1:length(names.true)){
  de <- possDe(t(adj.true), i, type = "dag")
  non.ancestor <- setdiff(1:length(names.true), de)
  anxx <- c(anxx, rep(i, length(de)))
  anyy <- c(anyy, de)
  xx <- c(xx, non.ancestor)
  yy <- c(yy, rep(i, length(non.ancestor)))
}
chosenXX.start <- names.true[anxx]
chosenYY.start <- names.true[anyy]

tmpstore <- matrix(0,100,2)

bg <- 1
  
chosenEdge <- edgeList[sample(nrow(edgeList), floor(nrow(edgeList)*bg)), ] 
positions <- sample(length(chosenXX.start), floor(length(chosenXX.start)*bg))
final.chosenXX <- chosenXX.start[positions]
final.chosenYY <- chosenYY.start[positions]

dc.pdag <- addBgKnowledge(res.cpdag, x = chosenEdge[, 1], y = chosenEdge[, 2])
res.dc.pdag.adj <- as(dc.pdag, 'matrix')

na.pdag <- add.bg.withlabel(res.cpdag, x = final.chosenXX, y = final.chosenYY)
res.na.pdag.adj <- as(na.pdag, 'matrix')

print(nrow(which(res.dc.pdag.adj == 1 & t(res.dc.pdag.adj) == 1, arr.ind = T))/2)
print(nrow(which(res.na.pdag.adj == 1 & t(res.na.pdag.adj) == 1, arr.ind = T))/2)

# ida
result.matrix <- matrix(0, nrow(res.dc.pdag.adj)*(nrow(res.dc.pdag.adj) -1 ), 15)
multiset.ida <- lapply(1:nrow(res.dc.pdag.adj), function(.) vector("list", nrow(res.dc.pdag.adj)))
multiset.dida <- lapply(1:nrow(res.dc.pdag.adj), function(.) vector("list", nrow(res.dc.pdag.adj)))
flag <- 1
for (i in 1:nrow(res.dc.pdag.adj)){
  for (j in 1:nrow(res.dc.pdag.adj)){
    if (i != j){
      cat('estimating casual effects for', i, 'and', j, '\n')
      # ida
      time.ida <- system.time(res.ida <- ida(i, j, cov(subsubdata), res.cpdag, method = "local", type = c("cpdag")))
      
      # sida
      time.sida <- system.time(res.sida <- ida(i, j, cov(subsubdata), na.pdag, method = "local", type = c("pdag")))
      
      # dida
      time.dida <- system.time(res.dida <- dida(i, j, cov(subsubdata), na.pdag))
      
      result.matrix[flag, ] <- c(time.ida[3], time.sida[3], time.dida[3], length(res.ida), length(res.sida), length(res.dida),
                                 min(res.ida), min(res.sida), min(res.dida), max(res.ida), max(res.sida), max(res.dida),
                                 mean(res.ida), mean(res.sida), mean(res.dida))
      multiset.ida[[i]][[j]] <- res.ida
      multiset.dida[[i]][[j]] <- res.dida
      cat('time elapsed',c(time.ida[3], time.sida[3], time.dida[3]), '\n')
      
      flag <- flag + 1
    }
  }
}

total.ida <- c()
total.dida <- c()
for (i in 1:nrow(res.dc.pdag.adj)){
  for (j in 1:nrow(res.dc.pdag.adj)){
    if (i != j){
      
      total.ida <- c(total.ida, multiset.ida[[i]][[j]])
      total.dida <- c(total.dida, multiset.dida[[i]][[j]])

      flag <- flag + 1
    }
  }
}


  # # ida
  # time.ida.dc <- system.time(res.ida.dc <- ida(pair[1], pair[2], cov(data), cpdag, method = "local", type = c("cpdag")))
  # 
  # # sida
  # time.sida.dc <- system.time(res.sida.dc <- ida(pair[1], pair[2], cov(data), pdag, method = "local", type = c("pdag")))
  # 
  # # dida
  # time.dida.dc <- system.time(res.dida.dc <- dida(pair[1], pair[2], cov(data), pdag))
  
  


# getting non-ancestral information. 
# for consistency, we only consider the ancestral relations implied by true dag

# consistency check
# final.chosenXX <- chosenXX.start
# final.chosenYY <- chosenYY.start
# for (i in 1:length(chosenXX)){
#   print(i)
#   if (sum(chosenXX.start == chosenXX[i] & chosenYY.start == chosenYY[i]) == 0){
#     # the considered edge is not contained in chosenXX.start
#     tmpx <- c(final.chosenXX, chosenXX[i])
#     tmpy <- c(final.chosenYY, chosenYY[i])
#     tmp.pdag <- add.bg.withlabel(res.cpdag, x = tmpx, y = tmpy)
#     if (is.null(tmp.pdag)){
#       cat(chosenXX[i], '--->', chosenYY[i], 'is not consistent. \n')
#     }else{
#       cat(chosenXX[i], '--->', chosenYY[i], 'is consistent. \n')
#       final.chosenXX <- tmpx
#       final.chosenYY <- tmpy
#     }
#   }
# }




















