
rm(list=ls())

library(igraph)
library(graph)
library(pcalg)
library(bnlearn)
library(CAM)
library(DREAM4)

pth <- getwd()
setwd(pth)
source(paste(pth, 'bgIDA.R', sep = "/", collapse = NULL))

load(paste(pth, 'result_direct_causal_real.RData', sep = "/", collapse = NULL))
load(paste(pth, 'result_non-ancestral_causal_real.RData', sep = "/", collapse = NULL))

set.seed(10)

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

final.ida <- matrix(0, length(r.ida), 4)
final.sida <- matrix(0, length(r.sida), 4)
final.dida <- matrix(0, length(r.dida), 4)

flag <- 1
for (nb in enb){
  
  for (bg in per){
    
    final.ida[flag, ] <- colMeans(r.ida[[flag]])
    final.sida[flag, ] <- colMeans(r.sida[[flag]])
    final.dida[flag, ] <- colMeans(r.dida[[flag]])
    
    flag <- flag + 1
  }

}


# plotting

# plot error bar
plot_error <- function(x, y, sd, len = 1, col = "black") { 
  len <- len * 0.05 
  arrows(x0 = x, y0 = y, x1 = x, y1 = y - sd, col = col, angle = 90, length = len) 
  arrows(x0 = x, y0 = y, x1 = x, y1 = y + sd, col = col, angle = 90, length = len) 
}

# ====================================== plot time ================================
c <- 1
flag <- 1
item <- c('time', 'average number of possible effects', 'average range of possible effects')

pdf(file = paste(c,'.pdf'), width = 10, height = 4)
par(mar = c(4,4,4,4))
par(oma = c(0.15,0.15,0.15,0.15))

x <- 1:length(final.ida[, 1])
x.up <- seq(0.5, length(final.ida[, 1])+0.5, length(per))
ylim <- c(0.09975, 2.11638)#range(final.ida[, c], final.sida[, c], final.dida[, c])*c(0.95, 1.05)
plot(x, final.ida[, 1], type = "n", xlab = "", ylab = "", ylim =ylim, xaxt = 'n')
# grid(col = "darkgrey", lty = 2)

#axis(1,at=x.up, label=rep('', 11))
axis(1,at=rep(1:11, 5) + 11*as.vector(t(matrix(rep(c(0, 2, 4, 6, 8), 11), 5, 11))),
     label=rep('', length(enb))[rep(1:11, 5) + 11*as.vector(t(matrix(rep(c(0, 2, 4, 6, 8), 11), 5, 11)))])
axis(1,at=rep(c(1,6,11), 5) + 11*as.vector(t(matrix(rep(c(0, 2, 4, 6, 8), 3), 5, 3))),
     label=rep(per, length(enb))[rep(c(1,6,11), 5) + 11*as.vector(t(matrix(rep(c(0, 2, 4, 6, 8), 3), 5, 3)))])

axis(3,at=rep(1:11, 5) + 11*as.vector(t(matrix(rep(c(1, 3, 5, 7, 9), 11), 5, 11))),
     label=rep('', length(enb))[rep(1:11, 5) + 11*as.vector(t(matrix(rep(c(1, 3, 5, 7, 9), 11), 5, 11)))])
axis(3,at=rep(c(1,6,11), 5) + 11*as.vector(t(matrix(rep(c(1, 3, 5, 7, 9), 3), 5, 3))),
     label=rep(per, length(enb))[rep(c(1,6,11), 5) + 11*as.vector(t(matrix(rep(c(1, 3, 5, 7, 9), 3), 5, 3)))])


axis(3,at=c(6, 17, 28, 39, 50, 61, 72, 83, 94, 105)[c(1,3,5,7,9)],
     label=c(expression(bold('en = 1')), expression(bold('en = 3')), expression(bold('en = 5')), expression(bold('en = 7')), expression(bold('en = 9'))))

axis(1,at=c(6, 17, 28, 39, 50, 61, 72, 83, 94, 105)[c(2,4,6,8,10)],
     label=c(expression(bold('en = 2')), expression(bold('en = 4')), expression(bold('en = 6')), expression(bold('en = 8')), expression(bold('en = 10'))))

polygon(c(x.up[1], x.up[2], x.up[2], x.up[1]), c(ylim[1], ylim[1], ylim[2], ylim[2]), col = "gray93", border = "black")
polygon(c(x.up[2], x.up[3], x.up[3], x.up[2]), c(ylim[1], ylim[1], ylim[2], ylim[2]), col = "white", border = "black")
polygon(c(x.up[3], x.up[4], x.up[4], x.up[3]), c(ylim[1], ylim[1], ylim[2], ylim[2]), col = "gray93", border = "black")
polygon(c(x.up[4], x.up[5], x.up[5], x.up[4]), c(ylim[1], ylim[1], ylim[2], ylim[2]), col = "white", border = "black")
polygon(c(x.up[5], x.up[6], x.up[6], x.up[5]), c(ylim[1], ylim[1], ylim[2], ylim[2]), col = "gray93", border = "black")
polygon(c(x.up[6], x.up[7], x.up[7], x.up[6]), c(ylim[1], ylim[1], ylim[2], ylim[2]), col = "white", border = "black")
polygon(c(x.up[7], x.up[8], x.up[8], x.up[7]), c(ylim[1], ylim[1], ylim[2], ylim[2]), col = "gray93", border = "black")
polygon(c(x.up[8], x.up[9], x.up[9], x.up[8]), c(ylim[1], ylim[1], ylim[2], ylim[2]), col = "white", border = "black")
polygon(c(x.up[9], x.up[10], x.up[10], x.up[9]), c(ylim[1], ylim[1], ylim[2], ylim[2]), col = "gray93", border = "black")
polygon(c(x.up[10], x.up[11], x.up[11], x.up[10]), c(ylim[1], ylim[1], ylim[2], ylim[2]), col = "white", border = "black")

i <- 0
for (nb in enb){
  
  x.sub <- (1:11) + i*11
  y1 <- matrix(0, 1, length(per))
  y2 <- matrix(0, 1, length(per))
  y3 <- matrix(0, 1, length(per))
  s1 <- matrix(0, 1, length(per))
  s2 <- matrix(0, 1, length(per))
  s3 <- matrix(0, 1, length(per))
  
  innerflag <- 1
  
  for (bg in per){
    
    if (c == 1 | c == 2){
      
      y1[innerflag] <- colMeans(r.ida[[flag]])[c]
      y2[innerflag] <- colMeans(r.sida[[flag]])[c]
      y3[innerflag] <- colMeans(r.dida[[flag]])[c]
      
      s1[innerflag] <- colSds(r.ida[[flag]])[c]
      s2[innerflag] <- colSds(r.sida[[flag]])[c]
      s3[innerflag] <- colSds(r.dida[[flag]])[c]
      
    }else{
      
      y1[innerflag] <- colMeans(r.ida[[flag]][,4] - r.ida[[flag]][,3])
      y2[innerflag] <- colMeans(r.sida[[flag]][,4] - r.sida[[flag]][,3])
      y3[innerflag] <- colMeans(r.dida[[flag]][,4] - r.dida[[flag]][,3])
      
      s1[innerflag] <- colSds(r.ida[[flag]][,4] - r.ida[[flag]][,3])
      s2[innerflag] <- colSds(r.sida[[flag]][,4] - r.sida[[flag]][,3])
      s3[innerflag] <- colSds(r.dida[[flag]][,4] - r.dida[[flag]][,3])
      
    }
    
    innerflag <- innerflag+1
    flag <- flag + 1
  }

  
  # pch = c(21, 22, 23, 24, 25)[mm]
  lines(x.sub, y1, type = "o", lty = 1, lwd = 1.5, pch = '', col = 'darkseagreen3')
  #plot_error(x.sub, y1, sd = s1, col = "black")
  
  points(x.sub, y2, type = "o",  lty = 1, lwd = 1.5, pch = '', col = 'cornflowerblue')
  #plot_error(x.sub, y2, sd = s1, col = "black")
  
  lines(x.sub, y3, type = "o",   lty = 1, lwd = 1.5, pch = '', col = 'coral')
  #plot_error(x.sub, y3, sd = s1, col = "black")
  
  i <- i+1

}


title(xlab = "percentage of background knowledge", ylab = item[c])

legend(x = c(1), y = c(2.075), legend = c('IDA', 'semi-local IDA', 'DIDA'), 
       cex = 0.8, lwd = 1.5, lty = 1, col = c('darkseagreen3','cornflowerblue','coral'), bg = 'white')
dev.off()


# =============================== plot size ========================================
c <- 2
flag <- 1
item <- c('time', 'average number of possible effects', 'average range of possible effects')

pdf(file = paste(c,'.pdf'), width = 10, height = 4)
par(mar = c(4,4,4,4))
par(oma = c(0.15,0.15,0.15,0.15))

x <- 1:length(final.ida[, 1])
x.up <- seq(0.5, length(final.ida[, 1])+0.5, length(per))
ylim <- c(0.950, 1.95)#range(final.ida[, c], final.sida[, c], final.dida[, c])*c(0.95, 1.1)
plot(x, final.ida[, 1], type = "n", xlab = "", ylab = "", ylim =ylim, xaxt = 'n')
# grid(col = "darkgrey", lty = 2)

#axis(1,at=x.up, label=rep('', 11))
axis(1,at=rep(1:11, 5) + 11*as.vector(t(matrix(rep(c(0, 2, 4, 6, 8), 11), 5, 11))),
     label=rep('', length(enb))[rep(1:11, 5) + 11*as.vector(t(matrix(rep(c(0, 2, 4, 6, 8), 11), 5, 11)))])
axis(1,at=rep(c(1,6,11), 5) + 11*as.vector(t(matrix(rep(c(0, 2, 4, 6, 8), 3), 5, 3))),
     label=rep(per, length(enb))[rep(c(1,6,11), 5) + 11*as.vector(t(matrix(rep(c(0, 2, 4, 6, 8), 3), 5, 3)))])

axis(3,at=rep(1:11, 5) + 11*as.vector(t(matrix(rep(c(1, 3, 5, 7, 9), 11), 5, 11))),
     label=rep('', length(enb))[rep(1:11, 5) + 11*as.vector(t(matrix(rep(c(1, 3, 5, 7, 9), 11), 5, 11)))])
axis(3,at=rep(c(1,6,11), 5) + 11*as.vector(t(matrix(rep(c(1, 3, 5, 7, 9), 3), 5, 3))),
     label=rep(per, length(enb))[rep(c(1,6,11), 5) + 11*as.vector(t(matrix(rep(c(1, 3, 5, 7, 9), 3), 5, 3)))])


axis(3,at=c(6, 17, 28, 39, 50, 61, 72, 83, 94, 105)[c(1,3,5,7,9)],
     label=c(expression(bold('en = 1')), expression(bold('en = 3')), expression(bold('en = 5')), expression(bold('en = 7')), expression(bold('en = 9'))))

axis(1,at=c(6, 17, 28, 39, 50, 61, 72, 83, 94, 105)[c(2,4,6,8,10)],
     label=c(expression(bold('en = 2')), expression(bold('en = 4')), expression(bold('en = 6')), expression(bold('en = 8')), expression(bold('en = 10'))))

polygon(c(x.up[1], x.up[2], x.up[2], x.up[1]), c(ylim[1], ylim[1], ylim[2], ylim[2]), col = "gray93", border = "black")
polygon(c(x.up[2], x.up[3], x.up[3], x.up[2]), c(ylim[1], ylim[1], ylim[2], ylim[2]), col = "white", border = "black")
polygon(c(x.up[3], x.up[4], x.up[4], x.up[3]), c(ylim[1], ylim[1], ylim[2], ylim[2]), col = "gray93", border = "black")
polygon(c(x.up[4], x.up[5], x.up[5], x.up[4]), c(ylim[1], ylim[1], ylim[2], ylim[2]), col = "white", border = "black")
polygon(c(x.up[5], x.up[6], x.up[6], x.up[5]), c(ylim[1], ylim[1], ylim[2], ylim[2]), col = "gray93", border = "black")
polygon(c(x.up[6], x.up[7], x.up[7], x.up[6]), c(ylim[1], ylim[1], ylim[2], ylim[2]), col = "white", border = "black")
polygon(c(x.up[7], x.up[8], x.up[8], x.up[7]), c(ylim[1], ylim[1], ylim[2], ylim[2]), col = "gray93", border = "black")
polygon(c(x.up[8], x.up[9], x.up[9], x.up[8]), c(ylim[1], ylim[1], ylim[2], ylim[2]), col = "white", border = "black")
polygon(c(x.up[9], x.up[10], x.up[10], x.up[9]), c(ylim[1], ylim[1], ylim[2], ylim[2]), col = "gray93", border = "black")
polygon(c(x.up[10], x.up[11], x.up[11], x.up[10]), c(ylim[1], ylim[1], ylim[2], ylim[2]), col = "white", border = "black")

i <- 0
for (nb in enb){
  
  x.sub <- (1:11) + i*11
  y1 <- matrix(0, 1, length(per))
  y2 <- matrix(0, 1, length(per))
  y3 <- matrix(0, 1, length(per))
  s1 <- matrix(0, 1, length(per))
  s2 <- matrix(0, 1, length(per))
  s3 <- matrix(0, 1, length(per))
  
  innerflag <- 1
  
  for (bg in per){
    
    if (c == 1 | c == 2){
      
      y1[innerflag] <- colMeans(r.ida[[flag]])[c]
      y2[innerflag] <- colMeans(r.sida[[flag]])[c]
      y3[innerflag] <- colMeans(r.dida[[flag]])[c]
      
      s1[innerflag] <- colSds(r.ida[[flag]])[c]
      s2[innerflag] <- colSds(r.sida[[flag]])[c]
      s3[innerflag] <- colSds(r.dida[[flag]])[c]
      
    }else{
      
      y1[innerflag] <- colMeans(r.ida[[flag]][,4] - r.ida[[flag]][,3])
      y2[innerflag] <- colMeans(r.sida[[flag]][,4] - r.sida[[flag]][,3])
      y3[innerflag] <- colMeans(r.dida[[flag]][,4] - r.dida[[flag]][,3])
      
      s1[innerflag] <- colSds(r.ida[[flag]][,4] - r.ida[[flag]][,3])
      s2[innerflag] <- colSds(r.sida[[flag]][,4] - r.sida[[flag]][,3])
      s3[innerflag] <- colSds(r.dida[[flag]][,4] - r.dida[[flag]][,3])
      
    }
    
    innerflag <- innerflag+1
    flag <- flag + 1
  }
  
  
  # pch = c(21, 22, 23, 24, 25)[mm]
  lines(x.sub, y1, type = "o", col = "darkseagreen3", lty = 1, lwd = 1.5, pch = '')
  #plot_error(x.sub, y1, sd = s1, col = "black")
  
  # lines(x.sub, y2, type = "o", col = "black", lty = 2, lwd = 1.3, pch = '')
  # #plot_error(x.sub, y2, sd = s1, col = "black")
  
  lines(x.sub, y3, type = "o",  col = "coral", lty = 1, lwd = 1.5, pch = '')
  #plot_error(x.sub, y3, sd = s1, col = "black")
  
  i <- i+1
  
}


title(xlab = "percentage of background knowledge", ylab = item[c])

legend(x = c(1), y = c(ylim[2])-0.015, 
       legend = c('without background knowledge (IDA)', 'with background knowledge (DIDA)'), 
       cex = 0.8, lwd = 1.5, lty = 1, bg = 'white', col = c('darkseagreen3', 'coral'))
dev.off()


# =============================== plot RANGE ========================================
flag <- 1
ratio <- lapply(1:110, function(.) c())
mean.ratio <- c()
std.ratio <- c()
for (nb in enb){
  
  for (bg in per){
    
    posi <- r.ida[[flag]][, 2] > 1
    
    diff.ida <- r.ida[[flag]][posi, 4] - r.ida[[flag]][posi, 3]
    diff.dida <- r.dida[[flag]][posi, 4] - r.dida[[flag]][posi, 3]
    ratio[[flag]] <- diff.dida/diff.ida
    mean.ratio <- c(mean.ratio, mean(ratio[[flag]]))
    std.ratio <- c(std.ratio, sd(ratio[[flag]]))
    
    flag <- flag + 1
  }
  
}

c <- 3
flag <- 1
item <- c('time', 'average number of possible effects', 'average range of possible effects')

pdf(file = paste(c,'.pdf'), width = 12, height = 4.5)
par(mar = c(4,4,4,4))
par(oma = c(0.15,0.15,0.15,0.15))

x <- 1:length(mean.ratio)
x.up <- seq(0.5, length(mean.ratio)+0.5, length(per))
ylim <- c(0, 1.5)
plot(x, mean.ratio, type = "n", xlab = "", ylab = "", ylim =ylim, xaxt = 'n')
# grid(col = "darkgrey", lty = 2)

#axis(1,at=x.up, label=rep('', 11))
axis(1,at=rep(1:11, 5) + 11*as.vector(t(matrix(rep(c(0, 2, 4, 6, 8), 11), 5, 11))),
     label=rep('', length(enb))[rep(1:11, 5) + 11*as.vector(t(matrix(rep(c(0, 2, 4, 6, 8), 11), 5, 11)))])
axis(1,at=rep(c(1,6,11), 5) + 11*as.vector(t(matrix(rep(c(0, 2, 4, 6, 8), 3), 5, 3))),
     label=rep(per, length(enb))[rep(c(1,6,11), 5) + 11*as.vector(t(matrix(rep(c(0, 2, 4, 6, 8), 3), 5, 3)))])

axis(3,at=rep(1:11, 5) + 11*as.vector(t(matrix(rep(c(1, 3, 5, 7, 9), 11), 5, 11))),
     label=rep('', length(enb))[rep(1:11, 5) + 11*as.vector(t(matrix(rep(c(1, 3, 5, 7, 9), 11), 5, 11)))])
axis(3,at=rep(c(1,6,11), 5) + 11*as.vector(t(matrix(rep(c(1, 3, 5, 7, 9), 3), 5, 3))),
     label=rep(per, length(enb))[rep(c(1,6,11), 5) + 11*as.vector(t(matrix(rep(c(1, 3, 5, 7, 9), 3), 5, 3)))])



axis(3,at=c(6, 17, 28, 39, 50, 61, 72, 83, 94, 105)[c(1,3,5,7,9)],
     label=c(expression(bold('en = 1')), expression(bold('en = 3')), expression(bold('en = 5')), expression(bold('en = 7')), expression(bold('en = 9'))))

axis(1,at=c(6, 17, 28, 39, 50, 61, 72, 83, 94, 105)[c(2,4,6,8,10)],
     label=c(expression(bold('en = 2')), expression(bold('en = 4')), expression(bold('en = 6')), expression(bold('en = 8')), expression(bold('en = 10'))))

polygon(c(x.up[1], x.up[2], x.up[2], x.up[1]), c(ylim[1], ylim[1], ylim[2], ylim[2]), col = "gray93", border = "black")
polygon(c(x.up[2], x.up[3], x.up[3], x.up[2]), c(ylim[1], ylim[1], ylim[2], ylim[2]), col = "white", border = "black")
polygon(c(x.up[3], x.up[4], x.up[4], x.up[3]), c(ylim[1], ylim[1], ylim[2], ylim[2]), col = "gray93", border = "black")
polygon(c(x.up[4], x.up[5], x.up[5], x.up[4]), c(ylim[1], ylim[1], ylim[2], ylim[2]), col = "white", border = "black")
polygon(c(x.up[5], x.up[6], x.up[6], x.up[5]), c(ylim[1], ylim[1], ylim[2], ylim[2]), col = "gray93", border = "black")
polygon(c(x.up[6], x.up[7], x.up[7], x.up[6]), c(ylim[1], ylim[1], ylim[2], ylim[2]), col = "white", border = "black")
polygon(c(x.up[7], x.up[8], x.up[8], x.up[7]), c(ylim[1], ylim[1], ylim[2], ylim[2]), col = "gray93", border = "black")
polygon(c(x.up[8], x.up[9], x.up[9], x.up[8]), c(ylim[1], ylim[1], ylim[2], ylim[2]), col = "white", border = "black")
polygon(c(x.up[9], x.up[10], x.up[10], x.up[9]), c(ylim[1], ylim[1], ylim[2], ylim[2]), col = "gray93", border = "black")
polygon(c(x.up[10], x.up[11], x.up[11], x.up[10]), c(ylim[1], ylim[1], ylim[2], ylim[2]), col = "white", border = "black")

i <- 0
for (nb in enb){
  
  x.sub <- (1:11) + i*11
  y <- mean.ratio[x.sub]
  s <- std.ratio[x.sub]
  

  
  
  # pch = c(21, 22, 23, 24, 25)[mm]
  lines(x.sub, y, type = "o", col = "coral", lty = 1, lwd = 1.5, pch = '')
  # plot_error(x.sub, y, sd = s, col = "black")
  
  # lines(x.sub, y2, type = "o", col = "black", lty = 2, lwd = 1.3, pch = '')
  # #plot_error(x.sub, y2, sd = s1, col = "black")
  
  # lines(x.sub, y3, type = "o",  col = "black", lty = 1, lwd = 1.3, pch = '')
  #plot_error(x.sub, y3, sd = s1, col = "black")
  
  i <- i+1
  
}


title(xlab = "percentage of background knowledge", ylab = item[c])

legend(x = c(1.2), y = c(ylim[2])-0.01, 
       legend = c('with background knowledge (DIDA)'), 
       cex = 0.8, lwd = 1.5, lty = c(1), bg = 'white', col = 'coral')
dev.off()

# ====================================== plot time detail ================================
p1 <- 12
p2 <- 17


# get data
t.ida.1 <- r.ida[[p1]][, 1]
t.ida.2 <- r.ida[[p2]][, 1]
t.sida.1 <- r.sida[[p1]][, 1]
t.sida.2 <- r.sida[[p2]][, 1]
t.dida.1 <- r.dida[[p1]][, 1]
t.dida.2 <- r.dida[[p2]][, 1]

# density fig1
d1 <- density(t.ida.1, bw = 0.025)
d2 <- density(t.sida.1, bw = 0.025)
d3 <- density(t.dida.1, bw = 0.025)

d4 <- density(t.ida.2, bw = 0.025)
d5 <- density(t.sida.2, bw = 0.025)
d6 <- density(t.dida.2, bw = 0.025)

xrange <- range(d1$x, d2$x, d3$x, d4$x, d5$x, d6$x)
yrange <- range(d1$y, d2$y, d3$y, d4$y, d5$y, d6$y)
xrange[1] <- 0
xrange[2] <- 1
yrange[1] <- 0

pdf(file = paste('density12','.pdf'), width = 4, height = 2.75)
par(mar = c(4,4,4,4))
par(oma = c(0.15,0.15,0.15,0.15))


x <- seq(0, xrange[2], length.out =100)
y <- seq(0, yrange[2], length.out = 100)
plot(x, y, type = "n", xlab = "", ylab = "", ylim =yrange, xaxt = 'n', yaxt = 'n')
grid(col = "darkgrey", lty = 2)

#axis(1,at=x.up, label=rep('', 11))
axis(1,at=seq(0, 1, 0.2), label=seq(0, 1, 0.2))
axis(2,at=seq(0, 12, 4), label=seq(0, 12, 4))

lines(d1, type = "o", lty = 1, lwd = 1.5, pch = '', col = 'darkseagreen3')

lines(d2, type = "o",  lty = 1, lwd = 1.5, pch = '', col = 'cornflowerblue')

lines(d3, type = "o",   lty = 1, lwd = 1.5, pch = '', col = 'coral')

title(xlab = "time", ylab = 'density')

legend(x = 'topright', legend = c('IDA', 'semi-local IDA', 'DIDA'), 
       cex = 0.8, lwd = 1.5, lty = 1, col = c('darkseagreen3','cornflowerblue','coral'), bg = 'white')
dev.off()

pdf(file = paste('density17','.pdf'), width = 4, height = 2.75)
par(mar = c(4,4,4,4))
par(oma = c(0.15,0.15,0.15,0.15))


x <- seq(0, xrange[2], length.out =100)
y <- seq(0, yrange[2], length.out = 100)
plot(x, y, type = "n", xlab = "", ylab = "", ylim =yrange, xaxt = 'n', yaxt = 'n')
grid(col = "darkgrey", lty = 2)

#axis(1,at=x.up, label=rep('', 11))
axis(1,at=seq(0, 1, 0.2), label=seq(0, 1, 0.2))
axis(2,at=seq(0, 12, 4), label=seq(0, 12, 4))

lines(d4, type = "o", lty = 1, lwd = 1.5, pch = '', col = 'darkseagreen3')

lines(d6, type = "o",   lty = 1, lwd = 1.5, pch = '', col = 'coral')

lines(d5, type = "o",  lty = 1, lwd = 1.5, pch = '', col = 'cornflowerblue')

title(xlab = "time", ylab = 'density')

legend(x = 'topright', legend = c('IDA', 'semi-local IDA', 'DIDA'), 
       cex = 0.8, lwd = 1.5, lty = 1, col = c('darkseagreen3','cornflowerblue','coral'), bg = 'white')
dev.off()

# ====================================== plot time detail 2================================
p1 <- 78
p2 <- 83


# get data
t.ida.1 <- r.ida[[p1]][, 1]
t.ida.2 <- r.ida[[p2]][, 1]
t.sida.1 <- r.sida[[p1]][, 1]
t.sida.2 <- r.sida[[p2]][, 1]
t.dida.1 <- r.dida[[p1]][, 1]
t.dida.2 <- r.dida[[p2]][, 1]

# density fig1
d1 <- density(t.ida.1, bw = 0.075)
d2 <- density(t.sida.1, bw = 0.075)
d3 <- density(t.dida.1, bw = 0.075)

d4 <- density(t.ida.2, bw = 0.075)
d5 <- density(t.sida.2, bw = 0.075)
d6 <- density(t.dida.2, bw = 0.075)

xrange <- range(d1$x, d2$x, d3$x, d4$x, d5$x, d6$x)
yrange <- range(d1$y, d2$y, d3$y, d4$y, d5$y, d6$y)
xrange[1] <- 0
xrange[2] <- 6
yrange[1] <- 0

pdf(file = paste('density78','.pdf'), width = 4, height = 2.75)
par(mar = c(4,4,4,4))
par(oma = c(0.15,0.15,0.15,0.15))


x <- seq(0, xrange[2], length.out =100)
y <- seq(0, yrange[2], length.out = 100)
plot(x, y, type = "n", xlab = "", ylab = "", ylim =yrange, xaxt = 'n', yaxt = 'n')
grid(col = "darkgrey", lty = 2)

#axis(1,at=x.up, label=rep('', 11))
axis(1,at=seq(0, 6, 2), label=seq(0, 6, 2))
axis(2,at=seq(0, 3, 0.5), label=seq(0, 3, 0.5))

lines(d1, type = "o", lty = 1, lwd = 1.5, pch = '', col = 'darkseagreen3')

lines(d2, type = "o",  lty = 1, lwd = 1.5, pch = '', col = 'cornflowerblue')

lines(d3, type = "o",   lty = 1, lwd = 1.5, pch = '', col = 'coral')

title(xlab = "time", ylab = 'density')

legend(x = 'topright', legend = c('IDA', 'semi-local IDA', 'DIDA'), 
       cex = 0.8, lwd = 1.5, lty = 1, col = c('darkseagreen3','cornflowerblue','coral'), bg = 'white')
dev.off()

pdf(file = paste('density83','.pdf'), width = 4, height = 2.75)
par(mar = c(4,4,4,4))
par(oma = c(0.15,0.15,0.15,0.15))


x <- seq(0, xrange[2], length.out =100)
y <- seq(0, yrange[2], length.out = 100)
plot(x, y, type = "n", xlab = "", ylab = "", ylim =yrange, xaxt = 'n', yaxt = 'n')
grid(col = "darkgrey", lty = 2)

#axis(1,at=x.up, label=rep('', 11))
axis(1,at=seq(0, 6, 2), label=seq(0, 6, 2))
axis(2,at=seq(0, 3, 0.5), label=seq(0, 3, 0.5))

lines(d4, type = "o", lty = 1, lwd = 1.5, pch = '', col = 'darkseagreen3')

lines(d5, type = "o",  lty = 1, lwd = 1.5, pch = '', col = 'cornflowerblue')

lines(d6, type = "o",   lty = 1, lwd = 1.5, pch = '', col = 'coral')

title(xlab = "time", ylab = 'density')

legend(x = 'topright', legend = c('IDA', 'semi-local IDA', 'DIDA'), 
       cex = 0.8, lwd = 1.5, lty = 1, col = c('darkseagreen3','cornflowerblue','coral'), bg = 'white')
dev.off()

# ====================================== plot time real data ================================

# get data
t.ida.1 <- result.matrix[, 1]
t.sida.1 <- result.matrix[, 2]
t.dida.1 <- result.matrix[, 3]


# density fig1
d1 <- density(t.ida.1, bw = 0.1)
d2 <- density(t.sida.1, bw = 0.1)
d3 <- density(t.dida.1, bw = 0.1)


xrange <- range(d1$x, d2$x, d3$x)
yrange <- range(d1$y, d2$y, d3$y)
xrange[1] <- 0
xrange[2] <- 8
yrange[1] <- 0

pdf(file = paste('real','.pdf'), width = 5, height = 3)
par(mar = c(4,4,4,4))
par(oma = c(0.15,0.15,0.15,0.15))


x <- seq(0, xrange[2], length.out =100)
y <- seq(0, yrange[2], length.out = 100)
plot(x, y, type = "n", xlab = "", ylab = "", ylim =yrange, xaxt = 'n', yaxt = 'n')
grid(col = "darkgrey", lty = 2)

#axis(1,at=x.up, label=rep('', 11))
axis(1,at=seq(0, 8, 2), label=seq(0, 8, 2))
axis(2,at=seq(0, 4, 1), label=seq(0, 4, 1))

lines(d1, type = "o", lty = 1, lwd = 1.5, pch = '', col = 'darkseagreen3')

lines(d2, type = "o",  lty = 1, lwd = 1.5, pch = '', col = 'cornflowerblue')

lines(d3, type = "o",   lty = 1, lwd = 1.5, pch = '', col = 'coral')

title(xlab = "time", ylab = 'density')

legend(x = 'topright', legend = c('IDA', 'semi-local IDA', 'DIDA'), 
       cex = 0.8, lwd = 1.5, lty = 1, col = c('darkseagreen3','cornflowerblue','coral'), bg = 'white')
dev.off()




  
