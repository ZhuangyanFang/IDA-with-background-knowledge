
rm(list=ls())

library(igraph)
library(graph)
library(pcalg)
library(bnlearn)
library(DREAM4)

pth <- getwd()
setwd(pth)
source(paste(pth, 'bgIDA.R', sep = "/", collapse = NULL))

set.seed(15)

# sample size
sampleSize <- 5000

# number of nodes
nn <- 50

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

for (nb in c(15)){
  
  for (bg in c(0)){
    
    r.ida[[length(r.ida) + 1]] <- matrix(0, sampletime, 4)
    r.sida[[length(r.sida) + 1]] <- matrix(0, sampletime, 4)
    r.dida[[length(r.dida) + 1]] <- matrix(0, sampletime, 4)
    
    for (st in 1:sampletime){
      cat('random graph, arbitrary , exp nbr size = ',nb, 'ratio of bg =',  bg, 'sample = ', st ,'\n')
      
      # sampling data
      dag <- pcalg::randomDAG(nn, nb/nn, lB = lb, uB = ub, V = as.character(1:nn))
      cpdag <- dag2essgraph(dag)
      data <- rmvDAG(sampleSize, dag, errDist="normal")
      
      # sampling non-ancestral background knowledge
      xx <- c()
      yy <- c()
      for (i in 1:nn){
        de <- possDe(t(as(dag, 'matrix')), i, type = "dag")
        non.ancestor <- setdiff(1:nn, de)
        xx <- c(xx, non.ancestor)
        yy <- c(yy, rep(i, length(non.ancestor)))
      }
      positions <- sample(length(xx), floor(length(xx)*bg))
      chosenXX <- xx[positions]
      chosenYY <- yy[positions]
      
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
      pdag <- add.bg(cpdag, x = chosenXX, y = chosenYY)
      # pdag <- as(pdag, 'graphNEL')
      
      
      save_dir = paste(pth, '\\samples','\\data','_',nn,'_',sampleSize,'_',nb,'_',bg,'_',st, '.RData', sep = "", collapse = NULL)
      save(list = c('dag','cpdag', 'data', 'chosenXX', 'chosenYY', 'pair'), file = save_dir)
      
      # ida
      time.ida <- system.time(res.ida <- ida(pair[1], pair[2], cov(data), cpdag, method = "local", type = c("cpdag")))
      
      # sida
      time.sida <- system.time(res.sida <- ida(pair[1], pair[2], cov(data), pdag, method = "local", type = c("pdag")))
      
      # dida
      time.dida <- system.time(res.dida <- dida(pair[1], pair[2], cov(data), pdag))
      
      if (sum(sort(unique(res.sida)) != sort(unique(res.dida)))!=0 ){
        break
      }
      
      r.ida[[length(r.ida)]][st,] <- c(time.ida[3], length(res.ida), min(res.ida), max(res.ida))
      r.sida[[length(r.sida)]][st,] <- c(time.sida[3], length(res.sida), min(res.sida), max(res.sida))
      r.dida[[length(r.dida)]][st,] <- c(time.dida[3], length(res.dida), min(res.dida), max(res.dida))
      
    } 
    
  }

}

save.image(paste(pth, 'result_non-ancestral_causal_real.RData', sep = "/", collapse = NULL))
