# Script to Generate Surrogate data simulating the calibration
# of a generic fluorescence probe in pollen tubes
# however the code attempts to mimick pHluorin data
# In: Portes, Damineli and Feijó (2021). Bio-protocol.
# "Spatiotemporal Quantification of Cytosolic pH in Arabidopsis Pollen Tubes"
# Code by Daniel Santa Cruz Damineli
# DISCLAIMER: CODE PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND.

##USER DEFINED SECTION##########################################################
# Define output directory
out.dir <- "~/Documents/bioprotocols/out/" # ATTENTION: Windows uses '\' instead of '/'

##AUTOMATIC SECTION: FUNCTION DEFINITIONS#######################################
################################################################################
GenerateSurrogateData <- function(out.dir = "~/Documents/bioprotocols/out/",
                                  pHs = rev(seq(5.5, 8, 0.25)),
                                  max.a2 = seq(20,85, length.out = length(pHs)),
                                  noise.ph = 0.25,
                                  rep.n = 9,
                                  channel.nms = c("FITC", "DAPI"),
                                  max.a1 = 20,
                                  tn = 50,
                                  X = seq(-20, 20, 0.1),
                                  y.intrcpt = 5,
                                  rate = 0.5,
                                  mu = 0, sig = 10000){
  
  # Define a linear model between pH and asymptote
  lmdl <- lm(max.a2 ~ pHs)
  
  # Initialize data frame to check data
  meds <- data.frame()
  
  # Loop thorugh pH values
  for(i in 1:length(pHs)){
    for(j in 1:rep.n){
      # Make a kymograph mimicking the unbound state wavelength
      kymo1 <- GenerateKymograph(tn = tn,
                                 X = X,
                                 max.a = max.a1,
                                 y.intrcpt = y.intrcpt,
                                 rate = rate, mu = mu, sig = sig)
      
      # Simulate experimental error in setting external pH
      actual.pH <- pHs[i] + rnorm(1, 0, noise.ph)
      actual.max <- lmdl$coefficients[2] * actual.pH + lmdl$coefficients[1]
      
      # Make a kymograph mimicking the bound state wavelength
      kymo2 <- GenerateKymograph(tn = tn,
                                 X = X,
                                 max.a = actual.max,
                                 y.intrcpt = y.intrcpt,
                                 rate = rate, mu = mu, sig = sig)
      
      # Store median ratios
      meds <- rbind(meds, 
                    cbind(pHs[i], j, 
                          median(kymo2[,10:110]/kymo1[,10:110], 
                                 na.rm = T)))
      
      # Generate file names following the convention specified in the manuscript
      fl.nm1 <- paste(channel.nms[1], "_", pHs[i], "_", j, ".txt", sep = "")
      fl.nm2 <- paste(channel.nms[2], "_", pHs[i], "_", j, ".txt", sep = "")
      
      # Save files
      colnames(kymo1) <- NULL
      write.table(kymo1, paste0(out.dir, fl.nm1), 
                  row.names = FALSE, col.names = FALSE)
      colnames(kymo2) <- NULL
      write.table(kymo2, paste0(out.dir, fl.nm2), 
                  row.names = FALSE, col.names = FALSE)
    }
  }
  
  return(meds)
}
################################################################################
GenerateKymograph <- function(tn = 50, 
                              X = seq(-20, 20, 0.1), 
                              max.a = 20, 
                              y.intrcpt = 5, 
                              rate = 0.5, mu = 0, sig = 1){
  # Generate an Asymptotic curve glued to a flat background
  sc <- GenerateMotherCurve(X, max.a, y.intrcpt, rate)
  # Make a Kymograph of the specified sized with growth of 1 pixel / frame
  kymo <- t(sapply(1:tn, 
                   function(i) sc[(1 + (tn - i)):(length(sc) - (i - 1))]))
  # Add noise
  noise <- t(sapply(1:tn, 
                    function(i) rnorm(length(sc) - tn + 1, mean = mu, sd = sig)))
  return(kymo + noise)
}
################################################################################
GenerateMotherCurve <- function(X = seq(-100, 20, 0.1), 
                                max.a = 20, 
                                y.intrcpt = 5, 
                                rate = 0.5){
  # Asymptotic curve
  Y <- max.a - (max.a - y.intrcpt) * exp (- rate * X)
  signal <- rev(Y) - min(Y)
  # Flat background
  bg <- rep(0, length(signal))
  return(c(signal, bg))  
}
################################################################################

##AUTOMATIC SECTION: CALL MAIN FUNCTION#########################################
# Run command and get sample output to see if it looks linear
if(!dir.exists(out.dir)){dir.create(out.dir)}
dat <- GenerateSurrogateData(out.dir = out.dir)

# Plot replicates and linear fit locally within R
plot(dat[, 1], dat[, 3], xlab = "pH", ylab = "DAPI/FITC")
abline(lm(dat[, 3] ~ dat[, 1]), col = "red")
