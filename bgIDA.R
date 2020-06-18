# bgIDA

library('pcalg')

criticalSet <- function(g, x, z){
  # given a chordal graph g, g may be disconnected
  # find the critical set of x w.r.t z
  # g is a graphNEL object
  # x is a number, z is an array
  
  n <- length(g@nodes)
  S <- list(c(x, 0, 0))
  C <- matrix(0, 1, n)
  flag <- 1
  dlp <- 0
  while (length(S) != 0){
    dlp <- dlp +1
    e <- S[[1]]
    S[[1]] <- NULL
    if (sum(z == e[1])==0 && e[2] != 0){
      for (alpha in setdiff(g@edgeL[[e[1]]][[1]],e[2])){
        if (sum(g@edgeL[[e[2]]][[1]]==alpha)==0){
          if (e[2] == x){
            S[[length(S)+1]] <- c(alpha, e[1], e[1])
          }else{
            S[[length(S)+1]] <- c(alpha, e[1], e[3])
          }
        }
      }
    }else if (sum(z == e[1])==0){
      for (alpha in g@edgeL[[e[1]]][[1]]){
        S[[length(S)+1]] <- c(alpha, e[1], alpha)
      }
    }else {
      C[flag] <- e[3]
      flag <- flag + 1
    }
    # if chordless cycle presents, then the while loop wil never end
    # DEAD loop prevent
    if (dlp > (n^2+1)) break
  }
  C <- unique(C[C!=0])
  return(C)
}

find.ancestors.and.chaincomp <- function(amat, x){
  # amat is a adj matrix of a cpdag
  # x is a given node
  # this function attempts to find the set of ancestors 
  # and chain components in the chain component
  
  n <- nrow(amat)
  amatSkel <- amat + t(amat)
  amatSkel[amatSkel != 0] <- 1
  amatDir <- amat - t(amat)
  amatDir[which(amatDir < 0)] <- 0
  amatUdir <- amat - amatDir
  
  w.an <- c(x)
  res.an <- c()
  while (length(w.an) > 0){
    node <- w.an[1]
    w.an <- w.an[-1]
    an.tmp <- which(amatDir[, node] == 1)
    w.an <- append(w.an, setdiff(setdiff(an.tmp, res.an), w.an))
    res.an <- append(res.an, node)
  }
  
  w.cp <- c(x)
  res.cp <- c()
  while (length(w.cp) > 0){
    node <- w.cp[1]
    w.cp <- w.cp[-1]
    cp.tmp <- which(amatUdir[, node] == 1)
    w.cp <- append(w.cp, setdiff(setdiff(cp.tmp, res.cp), w.cp))
    res.cp <- append(res.cp, node)
  }
  
  return(list(an = res.an, cp = res.cp))
  
}

add.bg <- function(cpdag, xx = c(), yy = c()){
  # adding non-ancestral background knowledge
  # since direct causal information is a special case of non-ancestral information
  # this function is a generalization of addBgKnowledge from pcalg
  # xx, yy are indices of nodes instead of labels
  
  amat <- as(cpdag, 'matrix')
  n <- nrow(amat)
  amatSkel <- amat + t(amat)
  amatSkel[amatSkel != 0] <- 1
  amatDir <- amat - t(amat)
  amatDir[which(amatDir < 0)] <- 0
  amatUdir <- amat - amatDir
  ug <- as(amatUdir, 'graphNEL')
  res.x <- c()
  res.y <- c()
  ancestor <- lapply(1:n, function(.) c())
  cp <- lapply(1:n, function(.) c())
  
  if (length(xx) != 0){
    for (i in 1:length(xx)){
      x <- xx[i]
      y <- yy[i]
      if (amatSkel[x, y] == 1){
        # x is adjacent to y
        res.x <- c(res.x, x)
        res.y <- c(res.y, y)
      }else{
        # x is not adjacent to y
        # y is not an ancestor of x
        # for each z \in an(x, cpdag) \cap cp(y), y is not an ancestor of z
        # find the critical set of y w.r.t z
        
        if (is.null(ancestor[[x]])){
          tmp <- find.ancestors.and.chaincomp(amat, x)
          ancestor[[x]] <- tmp$an
          cp[[x]] <- tmp$cp
        }
        if (is.null(cp[[y]])){
          tmp <- find.ancestors.and.chaincomp(amat, y)
          ancestor[[y]] <- tmp$an
          cp[[y]] <- tmp$cp
        }
        zset <- intersect(ancestor[[x]], cp[[y]])
        if (length(zset) != 0){
          c <- criticalSet(ug, y, zset)
          res.x <- c(res.x, c)
          res.y <- c(res.y, rep(y, length(c)))
        }
      }
    }
  }
  pdag <- addBgKnowledge(cpdag, x = res.x, y = res.y)
}

add.bg.withlabel <- function(cpdag, xx = c(), yy = c()){
  # another implementation of add.bg
  # if the nodes in cpdag have labels
  # and xx, yy are the labels of nodes instead of the indices of nodes
  # please use this function instead of add.bg
  
  amat <- as(cpdag, 'matrix')
  n <- nrow(amat)
  labels <- colnames(amat)
  amatSkel <- amat + t(amat)
  amatSkel[amatSkel != 0] <- 1
  amatDir <- amat - t(amat)
  amatDir[which(amatDir < 0)] <- 0
  amatUdir <- amat - amatDir
  ug <- as(amatUdir, 'graphNEL')
  # res.x and res.y store characters, i.e. nodes labels
  res.x <- c()
  res.y <- c()
  ancestor <- lapply(1:n, function(.) c())
  cp <- lapply(1:n, function(.) c())
  
  if (length(xx) != 0){
    for (i in 1:length(xx)){
      # x and y are labels
      x <- as.character(xx[i])
      y <- as.character(yy[i])
      pos.x <- which(labels == x)
      pos.y <- which(labels == y)
      if (amatSkel[x, y] == 1){
        # x is adjacent to y
        res.x <- c(res.x, x)
        res.y <- c(res.y, y)
      }else{
        # x is not adjacent to y
        # y is not an ancestor of x
        # for each z \in an(x, cpdag) \cap cp(y), y is not an ancestor of z
        # find the critical set of y w.r.t z
        
        if (is.null(ancestor[[pos.x]])){
          tmp <- find.ancestors.and.chaincomp(amat, pos.x)
          ancestor[[pos.x]] <- tmp$an
          cp[[pos.x]] <- tmp$cp
        }
        if (is.null(cp[[pos.y]])){
          tmp <- find.ancestors.and.chaincomp(amat, pos.y)
          ancestor[[pos.y]] <- tmp$an
          cp[[pos.y]] <- tmp$cp
        }
        zset <- intersect(ancestor[[x]], cp[[y]])
        if (length(zset) != 0){
          c <- criticalSet(ug, pos.y, zset)
          res.x <- c(res.x, labels[c])
          res.y <- c(res.y, rep(y, length(c)))
        }
      }
    }
  }
  pdag <- addBgKnowledge(cpdag, x = res.x, y = res.y)
}

dida <- function(x, y, mcov, graphEst, verbose = FALSE) {
  
  type = 'pdag'
  
  # graphEst is a maximal PDAG
  
  method <- "local"
  y.notparent = FALSE
  
  # graphEst <- addBgKnowledge(cpdag, x = bn[, 1], y = bn[, 2])
  amat <- ad.g <- wgtMatrix(graphEst)
  amat[which(amat != 0)] <- 1
  if (!isValidGraph(amat = amat, type = type)) {
    message("The input graph is not a valid ", type, ". See function isValidGraph() for details.\n")
  }
  nl <- colnames(amat)
  stopifnot(!is.null(nl))
  amatSkel <- amat + t(amat)
  amatSkel[amatSkel != 0] <- 1
  if (method == "local") {
    wgt.est <- (ad.g != 0) # cpdag, t
    tmp <- wgt.est - t(wgt.est)
    tmp[which(tmp < 0)] <- 0
    wgt.unique <- tmp # directed subgraph, t
    pa1 <- which(wgt.unique[x, ] != 0) # definite pa
    if (y %in% pa1) {
      beta.hat <- 0
    }
    else {
      wgt.ambig <- wgt.est - wgt.unique  # undirected subgraph
      pa2 <- which(wgt.ambig[x, ] != 0)  # sib
      if (verbose) 
        cat("\n\nx=", x, "y=", y, "\npa1=", pa1, "\npa2=", 
            pa2, "\n")
      if (length(pa2) == 0) {
        # calculate causal effect
        beta.hat <- lm.cov(mcov, y, c(x, pa1))
        if (verbose) 
          cat("Fit - y:", y, "x:", c(x, pa1), "|b.hat=", 
              beta.hat, "\n")
      }
      else {
        beta.hat <- NA
        ii <- 0
        pa2.f <- pa2
        pa2.t <- NULL
        if (hasExtensionNew(amat, amatSkel, x, pa1, pa2.t, 
                            pa2.f)) {
          ii <- ii + 1
          beta.hat[ii] <- lm.cov(mcov, y, c(x, pa1))
          if (verbose) 
            cat("Fit - y:", y, "x:", c(x, pa1), "|b.hat=", 
                beta.hat[ii], "\n")
        }
        for (i2 in seq_along(pa2)) {
          pa2.f <- pa2[-i2]
          pa2.t <- pa2[i2]
          if (hasExtensionNew(amat, amatSkel, x, pa1, pa2.t, 
                              pa2.f)) {
            ii <- ii + 1
            if (y %in% pa2.t) {
              beta.hat[ii] <- 0
            }
            else {
              beta.hat[ii] <- lm.cov(mcov, y, c(x, pa1, 
                                                pa2[i2]))
              if (verbose) 
                cat("Fit - y:", y, "x:", c(x, pa1, pa2[i2]), 
                    "|b.hat=", beta.hat[ii], "\n")
            }
          }
        }
        if (length(pa2) > 1) 
          for (i in 2:length(pa2)) {
            pa.tmp <- combn(pa2, i, simplify = TRUE)
            for (j in seq_len(ncol(pa.tmp))) {
              pa2.t <- pa.tmp[, j]
              pa2.f <- setdiff(pa2, pa2.t)
              if (hasExtensionNew(amat, amatSkel, x, pa1, 
                                  pa2.t, pa2.f)) {
                ii <- ii + 1
                if (y %in% pa2.t) {
                  beta.hat[ii] <- 0
                }
                else {
                  beta.hat[ii] <- lm.cov(mcov, y, c(x, 
                                                    pa1, pa2.t))
                  if (verbose) 
                    cat("Fit - y:", y, "x:", c(x, pa1, 
                                               pa2.t), "|b.hat=", beta.hat[ii], 
                        "\n")
                }
              }
            }
          }
      }
    }
  }
  unname(beta.hat)
}

hasExtensionNew <- function(amat, amatSkel, x, pa1, pa2.t, pa2.f){
  # borrowed from R packages pcalg
  # see, https://github.com/cran/pcalg
  
  res <- !has.new.coll.or.cycle(amat, amatSkel, x, pa1, pa2.t, pa2.f)
  res
}

has.new.coll.or.cycle <- function(amat,amatSkel, x, pa1, pa2.t, pa2.f) {
  ## Check if undirected edges that are pointed to x create a new v-structure
  ## or directed triangle containing x
  ## Additionally check, if edges that are pointed away from x create
  ## new v-structure; i.e. x -> pa <- papa would be problematic
  ## pa1 are definit parents of x
  ## pa2 are undirected "parents" of x
  ## pa2.t are the nodes in pa2 that are directed towards pa2
  ## pa2.f are the nodes in pa2 that are directed away from pa2
  ## Value is TRUE, if a new collider or a directed triangle is introduced 
  
  res <- FALSE
  
  # check v-structure
  if (length(pa2.t) > 0 && !all(is.na(pa2.t))) {
    ## check whether all pa1 and all pa2.t are connected;
    ## if not, there is a new collider
    if (length(pa1) > 0 && !all(is.na(pa1))) {
      res <- min(amatSkel[pa1, pa2.t]) == 0 ## TRUE if new collider
    }
    ## in addition, all pa2.t have to be connected
    if (!res && length(pa2.t) > 1) {
      A2 <- amatSkel[pa2.t,pa2.t]
      diag(A2) <- 1
      res <- min(A2) == 0 ## TRUE if new collider
    }
  }
  if (!res && length(pa2.f) > 0 && !all(is.na(pa2.f))) {
    ## consider here only the DIRECTED Parents of pa2.f
    ## remove undirected edges
    A <- amat-t(amat)
    A[A < 0] <- 0
    ## find parents of pa2.f
    cA <- colSums(A[pa2.f,,drop = FALSE])
    papa <- setdiff(which(cA != 0), x)
    ## if any node in papa is not directly connected to x, there is a new
    ## collider
    if (length(papa) > 0)
      res <- min(amatSkel[x,papa]) == 0 ## TRUE if new collider
  }
  
  ## checking direcrted triangle containing X
  if (!res ) {
    ## consider here pa = pa1 U pa2.t, cd = pa2.f U cd(amat)
    ## cd should not point towards pa
    
    ## adding new orientations to directed subgraph
    A <- amat-t(amat)
    A[A < 0] <- 0
    nA <- A
    nA[x, pa2.t] <- 1
    nA[pa2.t, x] <- 0
    nA[x, pa2.f] <- 0
    nA[pa2.f, x] <- 1
    
    # check whether cd point towards pa in nA
    cd <- which(nA[, x] != 0)
    pa <- which(nA[x, ] != 0) 
    if (length(cd) > 0 && length(pa) > 0){
      subA <- nA[pa, cd]
      res <- max(subA) == 1 ## TRUE if exists cd->pa
    }
  }
  res
}

lm.cov <- function (C, y, x) {
  # borrowed from R packages pcalg
  # see, https://github.com/cran/pcalg
  
  solve(C[x, x], C[x, y, drop = FALSE])[1, ]
}
