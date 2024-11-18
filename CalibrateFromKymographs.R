# Script to Generate a Calibration Linear fit from Kymographs 
# In: Portes, Damineli and Feijó (2021). Bio-protocol.
# "Spatiotemporal Quantification of Cytosolic pH in Arabidopsis Pollen Tubes"
# Accomplishes the steps on Section C. and video Supplementary File X
# See comments with STEP C.# in the main function to cross-reference with text
# Code by Daniel Santa Cruz Damineli
# DISCLAIMER: CODE PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND.

##USER DEFINED SECTION##########################################################
## Refer to Section D.1. of the manuscript for instructions
# Directory with Kymograph files that MUST follow:
fl.dir <- "~/Documents/bioprotocols/data/" # ATTENTION: Windows uses '\' instead of '/'
# Directory where the Calibration plot should be saved
out.dir <- "~/Documents/bioprotocols/out/" # ATTENTION: Windows uses '\' instead of '/'
# Specify the channel that should be in the numerator
top.ch <- "DAPI" # MUST be EXACTLY as written in the file prefix

# Margin from the tip to sample background 
mar <- 25 # in pixels
# Margin from the tip where ROI should start
roi.mar <- 25 # in pixels
# ROI width
roi.wdth <- 10 # in pixels
# Fraction of fluorescence value range to define threshold
frac <- 2

##AUTOMATIC SECTION: FUNCTION DEFINITIONS#######################################
# Main function where 'STEP D.2.' refers to Sections in the main manuscript
CalibrateFromKymographs <- function(fl.dir = "~/Documents/bioprotocols/data/",
                                    out.dir = "~/Documents/bioprotocols/out/",
                                    top.ch = "DAPI",
                                    frac = 2,
                                    mar = 25,
                                    roi.mar = 25,
                                    roi.wdth = 10){
  # Get file list from specified directory
  fls <- list.files(fl.dir)
  # Create output directory if it doesn't exist
  if(!dir.exists(out.dir)){dir.create(out.dir)}
  
  # Initiate empty data frame and vector of completed files
  ratio.dat <- data.frame()
  done.fls <- c()
  
  # Loop through files
  for(i in 1:length(fls)){
    fl <- fls[i]
    if(!(fl %in% done.fls)){
      # Split file to find matches
      fl.pieces <- strsplit(fl, split = "_")[[1]]
      experiment <- paste(fl.pieces[-1], collapse = "_")
      # Find matching experiment for both channels
      fl.pair <- fls[grep(experiment, fls)]
      
      if(length(fl.pair) != 2){
        # Error if there are not exactly 2 files
        message(paste("File pair not found for:", 
                      paste(fl.pair, collapse = " & ")))
      } else{
        message(paste("Analyzing:", 
                      paste(fl.pair, collapse = " & ")))
        
        # Parse files to get information separated by "_"
        fl.pair.info <- t(sapply(fl.pair, ParseFile))
        ch.ind <- which(fl.pair.info[, 1] == top.ch)
        
        ## STEP D.2.i. Open file
        kymo1 <- read.table(paste0(fl.dir, fl.pair[-ch.ind]), header = TRUE) + 0
        kymo1 <- t(apply(kymo1, 1, as.numeric))
        kymo2 <- read.table(paste0(fl.dir, fl.pair[ch.ind]), header = TRUE)
        kymo2 <- t(apply(kymo2, 1, as.numeric))
        
        ## STEP D.2.ii. Find Tip
        thrsh <- median(apply(kymo1, 1, 
                              function(vec) min(vec, na.rm = T) + 
                                (diff(range(vec, na.rm = T)) / frac)),
                        na.rm = T)
        
        pos <- apply(kymo1, 1, function(v) max(which(v > thrsh)))
        
        # Plot a single time point of the Kymograph to evaluate the threshold
        pdf(paste0(out.dir, strsplit(fl.pair[-ch.ind], split = "\\.")[[1]][1], 
                   "_KymoTimePoint_", floor(dim(kymo1)[1] / 2), ".pdf"), 
            height = 6, width = 8)
        plot(kymo1[floor(dim(kymo1)[1] / 2),], 
             xlab = "Pixel", ylab = "Fluorescence", 
             main = fl.pair[-ch.ind])
        abline(h = thrsh, col = "red", lwd = 2)
        abline(v = pos[floor(dim(kymo1)[1] / 2)], col = "green", lwd = 2)
        legend("topright", 
               legend = c(paste("Threshold =", formatC(thrsh, format = "e", digits = 2)), 
                          paste("Tip location =", pos[floor(dim(kymo1)[1] / 2)])), 
               col = c("red", "green"), text.col = c("red", "green"), lwd = 2, bty = "n")
        dev.off()
        
        ## STEP D.2.iii. Estimate background
        k1.bg <- median(RemoveOutliers(unlist(lapply(pos, function(p) 
          kymo1[, (p + mar):dim(kymo1)[2]]))))
        k2.bg <- median(RemoveOutliers(unlist(lapply(pos,function(p) 
          kymo2[, (p + mar):dim(kymo2)[2]]))))
        
        ## STEP D.2.iv-vi Calculate ratio
        ratio.kymo <- (kymo2 - k2.bg) / (kymo1 - k1.bg)
        ratio.kymo[is.infinite(ratio.kymo)] <- NA
        ratio.roi <- sapply(pos, function(p) 
          median(ratio.kymo[, (p - roi.mar - roi.wdth):(p - roi.mar - 1)], 
                 na.rm = T))
        
        pdf(paste0(out.dir, strsplit(fl.pair[-ch.ind], split = "\\.")[[1]][1], 
                   "_RatioInROI_", floor(dim(kymo1)[1] / 2), ".pdf"),
            height = 6, width = 8)

        par(mfrow=c(1, 2))
        plot(ratio.roi, type = "o", ylab = "Median ratio in ROI", xlab = "Time index")
        hist(ratio.roi, main = "Median ratio in ROI", 
             xlab = paste(fl.pair.info[ch.ind, 1], "/", fl.pair.info[-ch.ind, 1], sep = ""))
        dev.off()
        
        # Store in data frame
        ratio.dat <- rbind(ratio.dat, 
                           c(as.numeric(fl.pair.info[1, -1]), 
                             median(ratio.roi, na.rm = TRUE)))
        
        # Count files as done
        done.fls <- c(done.fls, fl.pair)
      }
    }
  }
  
  # Format ratio data
  ratio.dat <- as.data.frame(t(apply(ratio.dat, 1, as.numeric)))
  colnames(ratio.dat) <- c("Concentration", "Replicate", "Ratio")
  
  ## STEP D.2.vii. Linear fit of ratio data
  lfit <- lm(Ratio ~ Concentration, data = ratio.dat)
  lsumm <- summary(lfit)
  
  # Predict linear model to plot 95% CI
  conc <- sort(unique(ratio.dat$Concentration))
  conc.step <- max(diff(conc))
  xout <- seq(min(conc) - (2 * conc.step), 
              max(conc) + (2 * conc.step), 
              min(diff(conc)) / 10)
  lpred <- cbind(xout, 
                 predict(lfit, 
                         newdata = data.frame(Concentration = xout), 
                         interval = 'confidence'))
  
  
  # Data summary
  n.tbl <- table(ratio.dat$Concentration)
  se.tbl <- sapply(conc, 
                   function(cn) sd(ratio.dat$Ratio[ratio.dat$Concentration == cn]) / 
                     sqrt(n.tbl[names(n.tbl) == cn]))
  mean.tbl <- sapply(conc, 
                     function(cn) mean(ratio.dat$Ratio[ratio.dat$Concentration == cn]))
  
  # Plot Calibration fit
  pdf(paste0(out.dir, "CalibrationPlot.pdf"), height = 6, width = 8)
  plot(ratio.dat$Concentration, ratio.dat$Ratio, col = "grey", 
       ylab = colnames(ratio.dat)[3], 
       xlab = colnames(ratio.dat)[1])
  
  axis(3, at = sort(unique(ratio.dat$Concentration)))
  
  abline(lfit, col = "black", lwd = 2)
  polygon(y = c(lpred[, 4], rev(lpred[, 3])), x = c(xout, rev(xout)))
  points(conc, mean.tbl, col = "red", pch = 19)
  up <- mean.tbl + se.tbl
  down <- mean.tbl - se.tbl
  invisible(lapply(1:length(up), function(i) segments(x0 = conc[i], x1 = conc[i],
                                                      y0 = down[i], y1 = up[i],col = "red",lwd=2)))
  
  legend("topright", legend = c(paste0("R-squared = ", round(lsumm$adj.r.squared, digits = 2)),
                                paste0("p-val = ", formatC(lsumm$coefficients[2,][4], format = "e", digits = 2)),
                                paste0("Slope = ", round(lsumm$coefficients[2,][1], digits = 4)),
                                paste0("Intercept = ", round(lsumm$coefficients[1,][1], digits = 4))))
  
  text(x = conc, y = rep(min(ratio.dat$Ratio), length(n.tbl)), labels = paste0("n=", n.tbl), col = "blue", cex = 0.75)
  
  dev.off()
}
################################################################################
ParseFile <- function(fl){
  clmns <- strsplit(fl, split = "_")[[1]]
  clmns[length(clmns)] <- strsplit(clmns[length(clmns)], split = "\\.")[[1]][1]
  return(clmns)
}
################################################################################
MadOutliers <- function(vec, mad.tol = 3.5){
  mi <- 0.6745 * (vec - median(vec, na.rm = TRUE)) / mad(vec, na.rm = TRUE)
  return(which(abs(mi) > mad.tol))
}
################################################################################
RemoveOutliers <- function(vec, mad.tol = 3.5){
  return(vec[-MadOutliers(vec, mad.tol)])
}
################################################################################

##AUTOMATIC SECTION: CALL MAIN FUNCTION#########################################
CalibrateFromKymographs(fl.dir, out.dir, top.ch, frac, mar, roi.mar, roi.wdth)