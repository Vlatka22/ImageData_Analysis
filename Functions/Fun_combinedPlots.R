
#---------Function 'plotfit_minArea'-----------
# Plot fitted curves in time and 'cell-empty' area
plotfit_minArea <- function(fitvar, n1, n2, y_lab, colour, lpos, 
                            minArea = min_area,
                            fsize = 1, lw = 2,
                            miny, maxy, if_plot_min = T,
                            if_legend = T){
  par(mgp = c(1.5, .5, 0), mar = c(2.5, 2.5, 1.5, 5), bg = 'white') 
  # y-lim
  y_only <- unlist(lapply(fitvar[n1:n2], '[[', 2))
  if(missing(miny)) miny <- min(y_only)
  if(missing(maxy)) maxy <- max(y_only)
  # Plot set
  plot(1, type = 'n',
       xlim = c(0,  binNum),  ylim = c(miny,  maxy), xaxt = 'n', 
       xlab = 'x (um)', ylab = y_lab, main = '', cex.lab = fsize, cex.axis = fsize,
       col.axis = colour(1))
  axis(1, at = seq(0, 2000, 500) / (binW * umPix), labels = seq(0, 2000, 500),
       cex.axis = fsize)
  
  # Plot box demarking the cell density minima
  if(if_plot_min == T){
    lim <- par("usr")
    rect(minArea[1], lim[3], minArea[2], lim[4], border = "gray", col = "gray")
  }
  # Plot 
  cols = colour(n2-n1+1);
  for (i in n1:n2){
    lines(fitvar[[i]]$x, fitvar[[i]]$y, lwd = lw, col = cols[i-n1+1])
    if(i==frameNum) {break}
  }
  # Legend 
  # warning: "Increase right margin..."
  if(if_legend==T){
    suppressWarnings( gradientLegend(c(n1, i), color = cols, pos = lpos,
                                     coords = T, pos.num = NULL, n.seg = 1, dec = 0, 
                                     cex = fsize*.7, border.col = 'black') )
  }
}

#---------Function 'add_var2'-----------
# Overlay fitted curves in time
addvar2 <- function(fitvar2, n1, n2, col2, lpos2, fsize = 1, lw = 2, 
                    miny, maxy){
  # y-lim
  y_only <- unlist(lapply(fitvar2[n1:n2], '[[', 2))
  if(missing(miny)) miny <- min(y_only)
  if(missing(maxy)) maxy <- max(y_only)
  
  par(new = TRUE)
  plot(fitvar2[[1]]$x, fitvar2[[1]]$y, type = 'n', 
       axes = FALSE, xlab = '', ylab = '', ylim = c(miny, maxy))
  axis(4, cex.axis = fsize)
  
  # Add 2nd variable
  cols_2 = col2(n2-n1+1); 
  for (i in n1:n2){
    lines(fitvar2[[i]]$x, fitvar2[[i]]$y, lwd = lw, col = cols_2[i-n1+1])
    if(i==frameNum) {break}
  }
  # Legend
  # warning: "Increase right margin..."
  suppressWarnings( gradientLegend(c(n1, i), color = cols_2, pos = lpos2,
                                   coords = T, pos.num = NULL, n.seg = 1, dec = 0, 
                                   cex = fsize*.7, border.col = 'black') )
}

#---------Function 'plot_mix'-----------
# Plot expression, cell distribution (fitvar2) and 'cell-empty' area
plot_mix <- function(fitvar2, n1, n2, if_min = T, minArea = min_area, lw = 2,
                     if_legend1 = T, maxy2){
  # Plot expression
  plotfit_minArea(fit_maxI, n1, n2, '', colPal1, minArea = min_area,
                  lpos = pos_l_gene, lw = lw, if_legend = if_legend1)
  # y-lim
  y_only <- unlist(lapply(fitvar2[n1:n2], '[[', 2))
  if(missing(maxy2)) maxy2 <- max(y_only)
  
  # Plot cell distribution
  addvar2(fitvar2, n1, n2, colPal2, lpos2 = pos_l_nuc, lw = lw, maxy = maxy2)
}

#---------Function 'plot_mix_sp'-----------
# 'plot_mix', but with smoothness as an argument for both variables
# the 'box' is nuclei min detected by the default span parameter
plot_mix_sp <- function(n1, n2, span = .75, lw = 2, if_legend1 = T, maxy = maxy2){
  fitvar1 <- lw(maxI, span)
  fitvar2 <- lw(nuclei, span)
  #
  plotfit_minArea(fitvar1, n1, n2, '', colPal1, lpos = pos_l_gene,
                  lw = lw, if_legend = if_legend1)
  # y-lim
  y_only <- unlist(lapply(fitvar2[n1:n2], '[[', 2))
  if(missing(maxy)) maxy2 <- max(y_only)
  
  # Plot cell distribution
  addvar2(fitvar2, n1, n2, colPal2, lpos2 = pos_l_nuc, lw = lw, maxy = maxy2)
  title(main = span)
}

