
rm(list=ls())

library(igraph)
library(graph)
library(pcalg)
library(bnlearn)

pth <- getwd()
setwd(pth)
source(paste(pth, 'bgIDA.R', sep = "/", collapse = NULL))

set.seed(50)

# sample size
sampleSize <- 1000

# number of nodes
nn <- 100

# expected neigbour size
enb <- seq(1,10,1)

# percentage of background knowledge
per <- seq(0,1,0.1)

sampletime <- 100

alpha <- 0.01

lb <- 0.5
ub <- 2

r.ida <- list()
r.sida <- list()
r.dida <- list()

for (nb in c(10)){
  
  for (bg in c(0)){
    
    r.ida[[length(r.ida) + 1]] <- matrix(0, sampletime, 4)
    r.sida[[length(r.sida) + 1]] <- matrix(0, sampletime, 4)
    r.dida[[length(r.dida) + 1]] <- matrix(0, sampletime, 4)
    
    for (st in 1:sampletime){
      cat('random graph, exp nbr size = ',nb, 'ratio of bg =',  bg, 'sample = ', st ,'\n')
      
      # sampling data
      dag <- pcalg::randomDAG(nn, nb/nn, lB = lb, uB = ub, V = as.character(1:nn))
      cpdag <- dag2essgraph(dag)
      data <- rmvDAG(sampleSize, dag, errDist="normal")
      
      # sampling background knowledge
      amat <- as(dag, 'matrix')
      edgeList <- which(amat != 0, arr.ind = T)
      chosenEdge <- edgeList[sample(nrow(edgeList), floor(nrow(edgeList)*bg)), ]
      
      # sampling treatment and target
      pair <- sample(nn, 2)
      
      # # learning
      # p <- ncol(data)
      # data.suffStat <- list(C = cor(as.matrix(data)), n = nrow(data))
      # cat('learning graph \n')
      # pc <- pcalg::pc(suffStat = data.suffStat, indepTest = gaussCItest, alpha = alpha, p = p, skel.method = 'original')
      # cpdag.pc <- pc@graph
      # cat('pc finished \n')
      # pcs <- pc(suffStat = data.suffStat, indepTest = gaussCItest, alpha = alpha, p = p, skel.method = 'stable')
      # cpdag.pcs <- pcs@graph
      # cat('stable pc finished \n')
      
      # adding background knowledge
      pdag <- addBgKnowledge(cpdag, x = chosenEdge[, 1], y = chosenEdge[, 2])
      # pdag <- as(pdag, 'graphNEL')
      
      ida
      time.ida <- system.time(res.ida <- ida(pair[1], pair[2], cov(data), cpdag, method = "local", type = c("cpdag")))

      # sida
      time.sida <- system.time(res.sida <- ida(pair[1], pair[2], cov(data), pdag, method = "local", type = c("pdag")))

      # dida
      time.dida <- system.time(res.dida <- dida(pair[1], pair[2], cov(data), pdag))
      
      r.ida[[length(r.ida)]][st,] <- c(time.ida[3], length(res.ida), min(res.ida), max(res.ida))
      r.sida[[length(r.sida)]][st,] <- c(time.sida[3], length(res.sida), min(res.sida), max(res.sida))
      r.dida[[length(r.dida)]][st,] <- c(time.dida[3], length(res.dida), min(res.dida), max(res.dida))
      
    } 
  }
}
save.image(paste(pth, 'result_direct_causal_real.RData', sep = "/", collapse = NULL))

