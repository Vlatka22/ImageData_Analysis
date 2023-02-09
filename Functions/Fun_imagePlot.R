# --- Far outliers
outl <- function(x) {quantile(x, .75, na.rm = T) + 3 * IQR(x, na.rm = T)}

# --- Far outlier %
nout <- function(m) {sum(m > outl(m) & !is.na(m)) / sum(!is.na(m))}


#--------Function 'plot0'-----------
# Set the plot, convert NAs or zeros
plot0 <- function(m, repl_na = F, repl_0 = F, 
                  n1 = 1, n2 = dim(m)[1], tseq = 60,
                  boxbg = 'black', bg = 'gray87',
                  fsize = 1){
  par(mgp = c(1.6, .5, 0), mar = c(2.5, 2.5, 1.5, 5), bg = bg)
  
  # Set the plot, plot the box background (colour behind NAs)
  tmp <- m; tmp[,] <- 1
  image(1:ncol(m), 1:nrow(m), t(tmp), col = boxbg, xaxt = 'n', yaxt = 'n',
        xlab = 'x (um)', ylab = 'time (min)', ylim = c(n1, n2),
        cex.lab = fsize)
  axis(1, at = seq(0, 2000, 500) / (binW * umPix), labels = seq(0, 2000, 500),
       cex.axis = fsize)
  axis(2, at = seq(0, 300, tseq) / dt, labels = seq(0, 300, tseq), cex.axis = fsize)
  
  # NAs / zeros conversion (NAs are not plotted)
  # Replace original NAs with 0
  if (repl_na) {m <- m %>% replace(., is.na(.), 0)}
  # Replace original 0 to NAs 
  if (repl_0) {m <- m %>% replace(., .==0, NA)}
  
  return(m)
}


#--------Function 'plotmat'-----------
# Plot matrix as it is
plotmat <- function(m, colMap, title, repl_na = F, repl_0 = F, 
                    n1 = 1, n2 = dim(m)[1], tseq = 60,
                    fsize = 1.4, bg = 'gray87', # bg = 'white' for export
                    pos_legend = pos_l, dp = 0, if_legend = T){
  # Set the plot, convert NAs or zeros
  m <- plot0(m, repl_na, repl_0, n1, n2, tseq, bg, boxbg = 'black', fsize)
  title(main = title)
  
  # Plot
  image(1:ncol(m), 1:nrow(m), t(m), col = colMap, add = T)
  
  # Legend
  if (if_legend==T){
    gradientLegend(c(min(m, na.rm = T), max(m, na.rm = T)), color = colMap, 
                   pos = pos_legend, coords = T, n.seg = 2, dec = dp)
  }
}

#--------Function 'plotmat_outl'-----------
# Use the threshold as maximum for better image resolution
plotmat_outl <- function(m, colMap, title, repl_na = F, repl_0 = F,
                         n1 = 1, n2 = dim(m)[1], tseq = 60,
                         fsize = 1.4, bg = 'gray87', # bg = 'white' for export
                         pos_legend = pos_l, dp = 0, if_legend = T){
  # Set the plot, convert NAs or zeros
  m <- plot0(m, repl_na, repl_0, n1, n2, tseq, bg, boxbg = 'black', fsize)
  title(main = paste0(title, ' ...max:', round(max(m, na.rm = T))))
  
  # Define outlier threshold
  farout <- quantile(m, .75, na.rm = T) + 3 * IQR(m, na.rm = T)
  # Overwrite values above thr. with threshold value
  m <- m %>% replace(., . > farout, farout)
  
  # Plot
  image(1:ncol(m), 1:nrow(m), t(m), col = colMap, add = T)
  
  # Legend
  if (if_legend==T){
    gradientLegend(c(min(m, na.rm = T), max(m, na.rm = T)), color = colMap, 
                   pos = pos_legend, coords = T, n.seg = 2, dec = dp)
  }
}
