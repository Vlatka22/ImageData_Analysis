# --- Function lw
# loess (local polynomial regression fit):
lw <- function(binvar, sp = .75, repl_0 = F, repl_na = F){
  # sp: degree of smoothing
  fit <- list()
  # local var:
  frameNum = dim(binvar)[1]
  binNum = dim(binvar)[2]
  #
  if (repl_na) {binvar <- binvar %>% replace(., is.na(.), 0)}
  if (repl_0) {binvar <- binvar %>% replace(., .==0, NA)}
  #
  x = 1:binNum
  for (i in 1:frameNum){
    y = binvar[i,]
    tmp <- loess(y ~ x, span = sp)
    fit[[i]] <- list()
    fit[[i]]$x <- tmp$x #new x: e.g. places where y=NA removed etc.
    fit[[i]]$y <- tmp$fitted 
  } 
  return(fit)
}

#---------Plot dots & fitted curve-------------
plotfit <- function(binvar, fitvar, i, miny, maxy){
  # local var:
  binNum = dim(binvar)[2]
  #
  x = 1:binNum; y = binvar[i,] 
  nx = fitvar[[i]]$x; ny = fitvar[[i]]$y # btw. nx == x[!is.na(y)]
  # ylim
  if(missing(miny)) miny = min(binvar, na.rm = T)
  if(missing(maxy)) maxy = outl(binvar) #maxy = max(binvar, na.rm = T)
  # plot
  par(mgp = c(1.5, .5, 0), mar = c(2.5, 2.5, 1.5, 1.5), bg = 'gray87')
  plot(x, y, pch = 19, cex = .4, xaxt = 'n', xlab = 'x (um)',
       xlim = c(0, binNum),  ylim = c(miny, maxy))
  axis(1, at = seq(0, 2000, 500) / (binW * umPix), labels = seq(0, 2000, 500))
  lines(nx, ny, col = 'red', lwd = 3)
}

#---------Plot fitted curves in time-----------
plotfit_time <- function(fitvar, n1, n2, y_lab, miny, maxy,
                         colPal = viridis, if_legend = T){
  par(mgp = c(1.5, .5, 0), mar = c(2.5, 2.5, 1.5, 5))
  
  # y-lim
  y_only <- unlist(lapply(fitvar, '[[', 2))
  if(missing(miny)) miny <- min(y_only)
  if(missing(maxy)) maxy <- max(y_only)
  
  #Plot set
  plot(1, type = 'n',  bty='l',
       xlim = c(0,  binNum), ylim = c(miny,  maxy), xaxt = 'n', 
       xlab = 'x (um)', ylab = y_lab, main = '')
  axis(1, at = seq(0, 2000, 500) / (binW * umPix), labels = seq(0, 2000, 500))
  
  # Plot lines
  cols = colPal(n2-n1+1);
  for (i in n1:n2){
    lines(fitvar[[i]]$x, fitvar[[i]]$y, lwd = 2, col = cols[i-n1+1])
    if(i==frameNum) {break}
  }
  # Legend
  # warning: "Increase right margin..."
  suppressWarnings( gradientLegend(c(n1, i), color = cols, pos = .75, side = 4, dec = 0,
                 length = 0.25, depth = 0.02, inside = F, n.seg = 1) )
}

# ---  Compare different smoothness
sp_plot <- function(fit_a, fit_b, i, maxy){
  plot(fit_a[[i]]$y, pch = 19, cex = .5, ylim = c(0, maxy))
  points(fit_b[[i]]$y, pch = 19, cex = .5, col = 'purple')
  title(main = i)
}


