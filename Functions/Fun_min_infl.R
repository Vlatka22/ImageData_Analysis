
#---------Find curve's minimum-------------
find_min <- function(fitvar, i){
  x = fitvar[[i]]$x; y = fitvar[[i]]$y
  
  # y'= 0 & y''> 0
  # diff( diff(y) > 0) != 0 looks for a points of sign change (instead of y'=0,
  # as the data is discrete)
  tmp <- (diff( diff(y) > 0) != 0) & (diff(diff(y)/diff(x)) > 0)
  ext <- c(FALSE, tmp, FALSE) # append two els. 
  if(length(x[ext]) == 1) {return(x[ext])}
  else if (length(x[ext]) > 1) {
    ind = which(y[ext] == min(y[ext])) # if several minima, take the lowest
    return(x[ext][ind])
  } 
  else return(NA)
}

# --- Check detected minima
min_plot <- function(fit_n, min_n, i, lw = 10, maxy = maxy_n){
  par(mgp = c(1.4, .5, 0), mar = c(2.5, 2.5, 1.5, 1.5), bg = 'white')
  plot(fit_n[[i]]$y, pch = 19, cex = .75, ylim = c(0, maxy))
  abline(v = min_n[i], col = 'purple', lwd = lw)
}

# --- Function 'splitarea'
# Divide the area to: a) before the minimum and b) after the minimum
splitarea <- function(binvar, n1, n2, y_lab = '', y_lim, tseq = 60,
                      minArea = min_area, repl_NA = T,
                      col1 = 'lightskyblue3', col2 = 'purple', 
                      fsize = 1.4, lw = 4){
    # Replace NAs in binvar
  if (repl_NA == T){ 
    binvar <- binvar %>% replace(., is.na(.), 0)
  }
  
  # Mean values
  mean1 <- rowMeans(binvar[n1:n2, 1:minArea[1]], na.rm = T)
  mean2 <- rowMeans(binvar[n1:n2, minArea[2]:binNum], na.rm = T)
  
  # Plot
  par(mgp = c(1.5, .5, 0), mar = c(2.5, 2.5, 1.5, 2.5))
  plot(n1:n2, mean1, ylim = y_lim, type = 'l', lwd = lw,  col = col1,
       xaxt = 'n', xlab = 'time (min)', ylab = y_lab, cex.axis = fsize, cex.lab = fsize)
  axis(1, at = seq(0, 300, tseq) / dt, labels = seq(0, 300, tseq), cex.axis = fsize)
  lines(n1:n2, mean2, lwd = lw, col = col2)
}

#---------Find curve's inflection-----------
# In a specified direction (before or after the curve's maximum)
# Also, plot with detected extremes and biggest change, for control 
find_infl <- function(fitvar, i, win = 6, dd_sign = 'pos',
                      if_plot = T, if_result = F, 
                      csize = 1.5, lwidth = 3, y_lim = c(0, maxy_i)){
  # Simplify code
  x = fitvar[[i]]$x; y = fitvar[[i]]$y
  
  # Find maxima: y'=0 & y''<0
  tmp <- (diff( diff(y) > 0) != 0) & (diff(diff(y)/diff(x)) < 0)
  ext <- c(F, tmp, F) # add back two els
  
  # (dd_sign = 'pos'/'neg -> before/after highest maximum)
  max <- which(y == max(y[ext]))
  
  # Smooth the curve until there is only one inflection point
  ys = y; xs = x
  infl = c(T, T); ns = 0 # initiate
  while (sum(infl) > 1){
      ns = ns + 1
      ys <- rollapply(ys, win, mean)
      xs <- rollapply(xs, win, mean)
      # Inflection: y''==0
      # diff( diff(dy/dx) > 0) != 0 looks for a points of sign change of a first 
      # derivation (instead of y''=0, as the data is discrete)
      tmp <- diff( diff(diff(ys)/diff(xs)) > 0 ) != 0
      # add 2 els. & subset for increase/decrease
      # add 1 el. & subset for orientation to max
      if (dd_sign == 'pos'){
        infl <- c(F, F, tmp) & diff(ys) > 0
        if (length(max)==0) {infl <- c(F, infl)}
        else{
        infl <- c(F, infl) & xs < x[max]}
      } 
      else {
        infl <- c(F, F, tmp) & diff(ys) < 0
        if (length(max)==0) {infl <- c(F, infl)}
        else{
        infl <- c(F, infl) & xs > x[max]}
      }
  }
  # print(c(i,which(infl),max))
  
  # Plot
  if (if_plot){
    par(mgp = c(1.4, .5, 0), mar = c(2.5, 2.5, 1.5, 1.5))
    plot(x, y, pch = 19, cex = csize/2, col = 'darkgray', ylim = y_lim, 
         xlab = 'x', ylab = 'I', main = i)
    # Maxima
    points(x[ext], y[ext], pch = 19, cex = csize, col = 'blue')
    # Overlay smoothed curve
    lines(xs, ys, lwd = lwidth)
    # Inflection
    points(xs[infl] , ys[infl], pch = 19, cex = csize, col = 'purple') 
  }
  
  if (sum(infl)==0) {der = NA} else {der = xs[infl]}
  if(if_result) {return(der)}
}


#---------Find curve's inflection-----------
# In a specified direction (before or after the curve's maximum)
# Also, plot with detected extremes and biggest change, for control 
find_infl_0 <- function(fitvar, i, win = 6, dd_sign = 'pos',
                      if_plot = T, if_result = F, 
                      csize = 1.5, lwidth = 3, y_lim = c(0, maxy_i)){
  # Simplify code
  x = fitvar[[i]]$x; y = fitvar[[i]]$y
  
  # Find maxima: y'=0 & y''<0
  tmp <- (diff( diff(y) > 0) != 0) & (diff(diff(y)/diff(x)) < 0)
  ext <- c(F, tmp, F) # add back two els
  
  # (dd_sign = 'pos'/'neg -> before/after highest maximum)
  max <- which(y == max(y[ext]))
  
  # Smooth the curve until there is only one inflection point
  ys = y; xs = x
  infl = c(T, T); ns = 0 # initiate
  while (sum(infl) > 1){
    ns = ns + 1
    ys <- rollapply(ys, win, mean)
    xs <- rollapply(xs, win, mean)
    # Inflection: y''==0
    # diff( diff(dy/dx) > 0) != 0 looks for a points of sign change of a first 
    # derivation (instead of y''=0, as the data is discrete)
    tmp <- diff( diff(diff(ys)/diff(xs)) > 0 ) != 0
    # add 2 els. & subset for increase/decrease
    # add 1 el. & subset for orientation to max
    if (dd_sign == 'pos'){
      infl <- c(F, F, tmp) & diff(ys) > 0
      if (length(max)==0) {infl <- c(F, infl)}
      else {
        infl <- c(F, infl) & xs < x[max]
      }
    } else {
      infl <- c(F, F, tmp) & (diff(ys) < 0) 
      if (length(max)==0) {infl <- c(F, infl)}
      else {
        infl <- c(F, infl) & xs > x[max]
      }
    }
  }
  # print(c(i,which(infl),max))
  
  # Plot
  if (if_plot){
    par(mgp = c(1.4, .5, 0), mar = c(2.5, 2.5, 1.5, 1.5))
    plot(x, y, pch = 19, cex = csize/2, col = 'darkgray', ylim = y_lim, 
         xlab = 'x', ylab = 'I', main = i)
    # Maxima
    points(x[ext], y[ext], pch = 19, cex = csize, col = 'blue')
    # Overlay smoothed curve
    lines(xs, ys, lwd = lwidth)
    # Inflection
    points(xs[infl] , ys[infl], pch = 19, cex = csize, col = 'purple') 
  }
  
  if (sum(infl)==0) {der = NA} else {der = xs[infl]}
  if(if_result) {return(der)}
}