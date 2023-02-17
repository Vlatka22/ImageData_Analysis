library(matrixStats)
library(signal)
library(lmtest)
library(Cairo)

# All variables here are averaged over the coordinates in range2 
# (cells 'decided' to develop)

#---------Find flamindo minima-------------
# if_diff is T when fvar=dflam and F when fvar=flam
minima <- function(fvar, n1, n2, x_range = range2, win, y_lim, if_diff = F, 
                   if_plot = T, noise = -5){
  # y-lim
  if(missing(y_lim)) y_lim = c(min(fvar, na.rm = T), max(fvar, na.rm = T))
  
  # Find minima in time series defined as mean(fvar) in range2 over time n1-n2
  # Moving average with 'win' window length (dflam has sharp peaks: win = 1)
  x = rollapply(1:n2, width = win, mean)
  y = rollapply(rowMeans(fvar[1:n2, range2], na.rm = T), width = win, mean)
  
  # Find minima: y'=0 & y''>0
  # (y'=0: when does the sign of y' change, i.e. closest to zero)
  # using just diff(y) when all we need is sign...
  tmp <- (diff( diff(y) > 0) != 0) & (diff(diff(y)/diff(x)) > 0)
  ext <- c(F, tmp, F) # add back two els
  if (if_diff){
    ext <- ext & (y < noise) # remove noise from 'dflam' (smoothing destroys the peaks)
  }
  
  # Plot
  par(mgp = c(1.4, .5, 0), mar = c(2.5, 2.5, 1.5, 1.5), bg = 'white')
  if (if_plot){
    plot(x, y, xlim = c(n1, n2), ylim = y_lim, type = 'l', lwd = 2)
    # Minima
    abline(v = x[ext], col = 'steelblue', lty = 2, lwd = 1.5)
  }
  return(x[ext])
}

#---------Plot flamindo & gene expression dynamic-------------
flam_gene <- function(n1, n2, x_range = range2, minf = min_flam, lw = 5, fsize = 1.2, 
                      tseq = 15, y_lim = y_flam, y2_lim = y_int){
  # Mean values
  mean1 <- rowMeans(flam[1:n2, x_range], na.rm = T)
  mean2 <- rowMeans(maxI[1:n2, x_range], na.rm = T)
  
  # Plot flamindo
  par(mgp = c(1.5, .5, 0), mar = c(2.5, 2.5, 1.5, 2.5))
  plot(1:n2, mean1, type = 'l', lwd = lw,  xaxt = 'n', ylim = y_lim, 
       xlab = 'time (min)', ylab = 'Flamindo', cex.axis = fsize, cex.lab = fsize)
  axis(1, at = seq(0, 300, tseq) / dt, labels = seq(0, 300, tseq), cex.axis = fsize)
  
  # Flamindo minima
  abline(v = minf, col = 'steelblue', lty = 3, lwd = lw*.7)
  
  # Gene colour
  col2 = colPal1(3)[2]
  
  # Gene expression
  par(new = T)
  plot(1:n2, mean2, type = 'l', lwd = lw, col = col2, ylim = y2_lim,
       axes = FALSE, xlab = '', ylab = '')
  axis(4, cex.axis = fsize, col = col2, col.axis = col2)
}

#---------s golay filtering-------------
# adapted function to change y lims on plot 
# code from Elizabeth Westbrook (Jonathan Chubb's lab, LMCB UCL)
sgolay_line <- function(lw = 5, fsize = 1.2, tseq = 15, clim = y_flam, glimt = y_int,
                        add_flam_line = F, n = 3, dcAMP = F, 
                        dGene = F, p = 2, box = T, bg = 15){
  
  n <<- n
  p <<- p
  # Mean values
  
  
  var1 <- flam %>% replace(., is.na(.), 0)
  var2 <- maxI %>% replace(., is.na(.), 0)
  
  if (dcAMP == T && dGene == T){  mean1 <- sgolayfilt((rowMeans(var1[1:n3, range2])),p=p, n = n)
  mean2 <- sgolayfilt((rowMeans(var2[1:n3, range2])),p=p, n = n)
  
  bg1<- sgolayfilt((rowMeans(var1[1:n3, range2])),p=p, n = bg)
  bg2<- sgolayfilt((rowMeans(var2[1:n3, range2])),p=p, n = bg)
  
  x2 = c(1:n3)
  
  mean1 <- mean1 - bg1
  mean2 <- mean2 - bg2
  
  mean1 <- mean1*(-1)
  
  mean1 <- diff(mean1)/diff(x2)
  mean2 <- diff(mean2)/diff(x2)
  
  gtitle = paste0('d(',gene,')')
  ctitle = 'd(cAMP)'
  
  x3 = c(1:(n3-1))} 
  
  else if (dcAMP == F && dGene == T){
    mean1 <- sgolayfilt((rowMeans(var1[1:(n3-1), range2])),p=p, n = n)
    mean2 <- sgolayfilt((rowMeans(var2[1:n3, range2])),p=p, n = n)
    
    bg1<- sgolayfilt((rowMeans(var1[1:(n3-1), range2])),p=p, n = bg)
    bg2<- sgolayfilt((rowMeans(var2[1:n3, range2])),p=p, n = bg)
    
    x2 = c(1:n3)
    
    mean1 <- mean1 - bg1
    mean2 <- mean2 - bg2
    
    mean1 <- mean1*(-1)
    
    mean2 <- diff(mean2)/diff(x2)
    
    gtitle = paste0('d(',gene,')')
    ctitle = 'cAMP'
    
    x3 = c(1:(n3-1))} 
  
  else if (dcAMP == F && dGene == F){
    
    mean1 <- sgolayfilt((rowMeans(var1[1:n3, range2])),p=p, n = n)
    mean2 <- sgolayfilt((rowMeans(var2[1:n3, range2])),p=p, n = n)
    
    bg1<- sgolayfilt((rowMeans(var1[1:n3, range2])),p=p, n = bg)
    bg2<- sgolayfilt((rowMeans(var2[1:n3, range2])),p=p, n = bg)
    
    x2 = c(1:n3)
    
    mean1 <- mean1 - bg1
    mean2 <- mean2 - bg2
    
    mean1 <- mean1*(-1)
    x3 = x2
    gtitle = gene
    ctitle = 'cAMP'
    
  }
  
  par(mgp = c(1.5, .5, 0), mar = c(2.5, 2.5, 1.5, 2.5))
  # Flamindo
  plot(x3, mean1, type = 'l', lwd = lw, xaxt = 'n', yaxt = 'n', ylim = clim, 
       xlab = '', ylab = '', cex.axis = fsize, cex.lab = fsize) #, col = '#ff7f00'
  axis(1, at = seq(0, 300, tseq) / dt, labels = seq(0, 300, tseq), cex.axis = fsize)
  # axis(2, cex.axis = fsize, col = colPal1(100)[50], col.axis = colPal1(100)[50])#
  axis(4, cex.axis = fsize, ylab = ctitle) #col.axis = '#ff7f00'
  # Gene expression
  par(new = T)
  plot(x3, mean2, type = 'l', lwd = lw, ylim = glim,
       col = colPal1(100)[50], axes = FALSE, xlab = 'time (min)', ylab = gtitle)
  axis(2, cex.axis = fsize, col =colPal1(100)[50] , col.axis = colPal1(100)[50]) # colPal1(100)[50]
  mtext(ctitle, 4, line =1.3)
  
  if(add_flam_line == T){
    abline(v = min_flam, col = 'steelblue', lty = 3, lwd = lw*.7)
  }
  #abline(v = min_flam, col = 'steelblue', lty = 3, lwd = lw*.7)
  
  if(box == T){
    box()
  }
}


#----------------------
# adapted function to change y lims on plot 
# code from Elizabeth Westbrook (Jonathan Chubb's lab, LMCB UCL)
plot_ccf <- function(n = 3, p = 1, bg = 15, lag = 5){
  # Mean values
  
  var1 <- flam %>% replace(., is.na(.), 0)
  var2 <- maxI %>% replace(., is.na(.), 0)
  
  
  mean1 <- sgolayfilt((rowMeans(var1[1:n2, range2])),p=p, n = n)
  mean2 <- sgolayfilt((rowMeans(var2[1:n2, range2])),p=p, n = n)
  
  bg1<- sgolayfilt((rowMeans(var1[1:n2, range2])),p=p, n = bg)
  bg2<- sgolayfilt((rowMeans(var2[1:n2, range2])),p=p, n = bg)
  
  x2 = c(1:n2)
  
  mean1 <- mean1 - bg1
  mean2 <- mean2 - bg2
  
  mean1 <- mean1*(-1)
  
  #mean2 <- diff(mean2)/diff(x2)
  
  x3 = c(1:(n2-1))
  
  par(mgp = c(1.5, .5, 0), mar = c(2.5, 2.5, .5, .5))
  return(list(ccf(mean1, mean2,lag)))
  
}
