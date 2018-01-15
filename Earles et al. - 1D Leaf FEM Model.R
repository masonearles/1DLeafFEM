###############################################################
# 1-D Steady-state Porous-media Diffusion Reaction Leaf Model #
# for incorporation of 1-D absorption profiles                #
# Author: J. Mason Earles et al.                              #
# Updated: Feb. 27th, 2017                                    #
###############################################################

# Load libraries
library(deSolve)
library(rootSolve)
library(ReacTran)
library(ggplot2)
library(gridExtra)
library(foreach)
library(doParallel)
library(Cairo)
library(cowplot)
library(colorspace)
library(caTools)
library(plyr)

# (If multicore processing is desired) Set up cluster
cl = makeCluster(detectCores())
registerDoParallel(cl)
#stopCluster(cl) # Use this command to stop cluster following simulation

# Load Brodersen et al. (2008) light response data
# Set your local working directory to '.../Script'
d.mu = read.table('_inputData/Brodersen2008SummaryData.txt')

# Load absorption profile data
dir_inputData = '_inputData/'
LeafOut.List = list()
LIPNames = dir(dir_inputData)[grep(pattern="LightIntensityProfile",dir(dir_inputData))]
for(n in 1:length(LIPNames)) {
alpha_table = read.table(paste0(dir_inputData,LIPNames[n]), header=TRUE) # Light absorption profile
alpha_table[,1] = rev(alpha_table[,1]); alpha_table[,2] = rev(alpha_table[,2]) #Reverse light absorption profile

# Define simulation parameters
nx = N = dim(alpha_table)[1] # Number of spatial steps
t_leaf = alpha_table[nx,1] # Leaf thickness; also mesophyll thickness [m]
t_elem = t_leaf/nx; # Element size [m]
t_elem_um = t_elem*1e6 # Element size [um]
Intdepth = seq(0,t_leaf, by = t_elem) # Positions of element interfaces [m]
Nint = length(Intdepth) # Number of element interfaces
Depth = seq(t_elem,t_leaf,by=t_elem) # Positions of element centroids [m]
moltoppm = 400/.0172 # mol to ppm conversion for CO2

# Generate  vector with multiple irradiance intensities 
I.vec <- c(seq(0, 0.5e-3, 5e-5),0.6e-3,0.75e-3,1e-3,1.25e-3,1.5e-3,2e-3,5e-3) # Irradiance levels of interest [mol m-2 s-1]
leaf.vec <- 1:length(I.vec) # Vector for storing output

# Define values for parameteric sensitivity analysis; Parameters defined below in comments
scenario.out <- list()
.g_liq = c(0.25e-3,1e-3,1e-2)
.frac_pal = if(grepl('Sun',LIPNames[n])) c(0.6,0.45,0.75) else c(0.45,0.3,0.6)
.por_spg = c(0.3,0.1,0.5)
.por_pal = c(0.1,0.05,0.15)
.tort = c(1.55,1.15,1.95)
.Sm_spg = c(5,2.5,7.5)
.Sm_pal = c(30,15,35)
.Vstrom = c(1.74e-6,1.24e-6,2.24e-6)
.Vmito = 0.27e-7
.Alpha = if(grepl('Sun',LIPNames[n])) c(0.72,0.4,1.0) else c(0.69,0.4,1.0)
.Alpha = if(grepl('Direct',LIPNames[n])) 1*.Alpha else 0.96*.Alpha
.Beta = c(0.44,0.4,0.5)
.k_c = c(3,2,4)
.X_c = if(grepl('Sun',LIPNames[n])) c(2.5,1.5,3.5) else c(1,0.5,1.5)
.K_m = c(18.7e-3,12.7e-3,24.7e-3)
.Gamma = 1.75e-3
.Theta = c(1,0.3,0.9);
.Rd = 0.066;
.Jmax = if(grepl('Sun',LIPNames[n])) c(275e-6,150e-6,350e-6) else c(275e-6,150e-6,350e-6)

# Define scenario for testing spatial distribution of phiPSII; Uncomment scenario of interest
phiPSII = if(grepl('Sun',LIPNames[n])) seq(0.85,0.5,length.out=nx) else seq(0.5,0.5,length.out=nx) # Testing assumption that quantum efficiency of PSII varies from 0.5 at the adaxial surface to 0.85 at the abaxial surface
#phiPSII = if(grepl('Sun',LIPNames[n])) seq(0.85,0.85,length.out=nx) else seq(0.85,0.85,length.out=nx) # Testing assumption that quantum efficiency of PSII constant at 0.85 throughout leaf
#phiPSII = if(grepl('Sun',LIPNames[n])) seq(0.5,0.5,length.out=nx) else seq(0.5,0.5,length.out=nx) # Testing assumption that quantum efficiency of PSII constant at 0.5 throughout leaf

# Create a dataframe with parameters of interest for sensitivity analysis
parameters.df = expand.grid(list(.g_liq = .g_liq[1], .frac_pal = .frac_pal[1], .por_spg = .por_spg[1], .por_pal = .por_pal[1], .tort = .tort[1], .Sm_spg = .Sm_spg[1], .Sm_pal = .Sm_pal[1], .Vstrom = .Vstrom[1], .Vmito = .Vmito[1], .Alpha = .Alpha[1], .Beta = .Beta[1], .k_c = .k_c[1], .X_c = .X_c[1], .K_m = .K_m[1], .Gamma = .Gamma[1], .Theta = .Theta[1], .Rd = .Rd[1], .Jmax = .Jmax[1]))
scenarios.df = parameters.df

# Add scenarios to baseline scenario dataframe 
scenarios = c('.g_liq','.frac_pal','.por_spg','.por_pal','.Sm_spg','.Sm_pal','.Vstrom','.tort','.Alpha','.Beta','.k_c','.X_c','.K_m','.Jmax','.Theta')
for(Y in scenarios){
  for(k in 2:3){
    parameters.df[,which(colnames(parameters.df)==Y)] = get(Y)[k]
    scenarios.df = rbind(scenarios.df,parameters.df)
    parameters.df[,which(colnames(parameters.df)==Y)] = get(Y)[1]
  }
}
scenarios.df$PARM = c('baseline',paste0(rep(scenarios,each=2))) #Column that defines which parameter was varied

#########################################################################
# Run 1-D reaction porous media reaction-diffusion photosynthesis model #
#########################################################################

scenario.out <- foreach(j=1:nrow(scenarios.df), .export='steady.1D') %dopar% {
  
  # Geometric variables
  frac_pal = scenarios.df$.frac_pal[j] # Fraction palisade mesophyll
  frac_spg = 1-frac_pal # Fraction spongy mesophyll
  Sm_spg = scenarios.df$.Sm_spg[j] # Sm spongy mesophyll [m2 m-2]
  Sm_pal = scenarios.df$.Sm_pal[j] # Sm palisade mesophyll [m2 m-2]
  Sm = c(rep(Sm_spg/(nx+1),floor((nx+1)*frac_spg)),rep(Sm_pal/(nx+1),ceiling(nx*frac_pal)))[1:nx]*(t_elem*1e6) # Surface area mesophyll per surface area leaf [m2 m-2]
  Vstrom = scenarios.df$.Vstrom[j] # Stroma volume per mesophyll surface area [m3 m-2]
  Vmito = scenarios.df$.Vmito[j] # Mitochondrial volume per mesophyll surface area [m3 m-2]
  V_mito = Vmito/Vstrom # Mitochondrial volume per stroma volume [m3 m-3]
  por_spg = scenarios.df$.por_spg[j] # Porosity of the spongy mesophyll [m3 air m3 leaf]
  por_pal = scenarios.df$.por_pal[j] # Porosity of the palisade mesophyll [m3 air m3 leaf]
  por = c(rep(por_spg,floor((nx+1)*frac_spg)),rep(por_pal,ceiling(nx*frac_pal)))[1:nx] # Porosity [m3 air m3 leaf]
  tort = scenarios.df$.tort[j] # Tortuosity of the palisade and spongy mesophyll [m m-1]
  g_liq = scenarios.df$.g_liq[j] # Cell wall + liquid conductivity into stroma [m s-1]
  Intpor = c(por,por[N]) # Porosity at interface [m3 m-3]
  Dz = 1.54e-5 # Diffusivity of CO2 in air [m2 s-1]
  
  # Photosynthetic variables
  moltoppm = 400/.0172
  Alpha = scenarios.df$.Alpha[j] # Leaf level absorption [mol mol-1]
  Beta = scenarios.df$.Beta[j] # Fraction of light absorbed by PSII [mol mol-1]
  k_c = scenarios.df$.k_c[j] # Rubisco turnover rate [s-1] 
  X_c = scenarios.df$.X_c[j] # Rubisco concentration [mol m-3]
  K_m = scenarios.df$.K_m[j] # Rubisco effective Michaelis-Menten constant [mol m-3]
  Gamma = scenarios.df$.Gamma[j] # CO2 compensation point [mol m-3]
  Theta = scenarios.df$.Theta[j] # Curvature factor in FvCB model
  Rd = scenarios.df$.Rd[j] # Dark respiratory rate [mol m-3 s-1]
  Sm_std = 30 # Sm at which assumed J_max occurs
  
  ##
  # Set Rubisco concentration based on Nishio et al. (1993)
  ##
  Xc_dist_z = seq(1,1,nx)
  Xc_dist_z = Xc_dist_z*c(seq(0.75,0.75,length.out=round(0.28*nx)),seq(0.75,1.35,length.out=round(0.43*nx)),seq(1.35,0.85,length.out=round(0.29*nx)))
  X_c_z = X_c * Xc_dist_z
  
  ##
  # Define scenario for testing spatial distribution of Jmax; Uncomment scenario of interest
  ##
  
  Jprop = scenarios.df$.Jmax[j]/Sm_std # Proportionality factor between mesophyll surface area and maximum e- transport rate
  J_max = Jprop*sum(Sm) # Maximum e- transfer rate on a leaf area basis [mol m-2 s-1]
  j_max_z = J_max/sum(Sm)/Vstrom # Electron transport rate on a stroma volume basis at position z [mol m-3 s-1]
  
  # Scenario 1: Jmax proportional to Rubisco distribution
  j_dist_z = c(seq(0.4,0.4,length.out=0.28*nx),seq(0.4,1,length.out=0.43*nx),seq(1,0.8,length.out=0.29*nx))[1:nx]
  j_max_z_2 = j_max_z*j_dist_z/sum(j_dist_z)
  
  # Scenario 2: Jmax proportional to chloroplast distribution
  #j_max_z_2 = j_max_z*alpha_table[,2]/sum(alpha_table[,2])
  
  # Scenario 3: Jmax proportional to stroma volume
  #j_max_z_2 = j_max_z*Sm/sum(Sm)
  
  # Scenario 4: Jmax proportional to chloroplast distribution
  #j_dist_z = c(seq(0,1,length.out=0.2*nx),seq(1,1,length.out=0.2*nx),seq(1,0.8,length.out=0.2*nx),seq(0.8,0.4,length.out=0.25*nx),seq(0.4,0,length.out=0.15*nx))[1:nx]
  #j_max_z_2 = j_max_z*j_dist_z/sum(j_dist_z)

  # Set scenario
  j_max_z = j_max_z_2
  
  ##
  # Generate light response curves for multiple irradiance intensities
  ##
  out.list <- list() # List for storing output
  An.out <- list() # List for storing output
  j.out <- list() # List for storing output

  for(i in 1:length(I.vec)){
    
    # Environmental variables
    I_0 = I.vec[i] # Incident irradiance [mol m-2 s-1]
    
    # Define electron transport rate function
    I_z = I_0*Beta*Alpha*phiPSII*(alpha_table[,2]/sum(alpha_table[,2]))/sum(Sm)/Vstrom
    j_z = pmin(I_z,j_max_z)

    # Set boundary conditions and initial conditions
    C_s = 0.85*1.72e-2 # CO2 concentration at stomate
    C_ias = rep(0.5*1.72e-2, nx) # Define initial condition of [CO2] in IAS
    C_liq = C_ias*0.5 # Define initial condition of [CO2] in chloroplast stroma
    y <- c(C_ias,C_liq) # Setup [CO2] vector for vapor and liquid
    
    
    ######################################
    # Define reaction-diffusion function #
    ######################################
    
    DifRxn <- function(t, y, parms) {
      C_ias <- y[1:N] # Define vapor [CO2]
      C_liq <- y[(N+1):(2*N)] # Define liquid [CO2]
      
      Flux_Cias <- -Dz*(Intpor/tort)*diff(c(C_s,C_ias,C_s))/t_elem # Calculate change in vapor [CO2] over time; set constant [CO2] = C_s at lower/upper boundary of leaf
      
      #Flux_Cias[N] <- 0 # Set no flux boundary condition at top of leaf; for hypostomatous leaf
      
      An <- pmin((k_c*X_c_z*C_liq)/(K_m+C_liq),nx*(C_liq*j_z/(4*C_liq+8*Gamma))) # Calculate carboxylation rate
      
      Rp <- (An*Gamma)/C_liq # Calculate oxygenation rate
      
      dC_ias <- (-diff(Flux_Cias)/t_elem + g_liq*(C_liq - C_ias)/t_elem)/por # Calculate the divergence of the C_ias flux
      dC_liq <- (g_liq*(C_ias - C_liq)/t_elem - An + Rp + Rd)*t_leaf/sum(Sm)/Vstrom # Caculate the corresponding change in C_liq flux
      
      return(list(c(dC_ias, dC_liq))) # return output from function
    }
    
    # Solve 1-D Steady-state Porous Media Leaf Reaction-Diffusion Model
    out <- steady.1D(y = y, func = DifRxn, parms = NULL, nspec = 2, names = c("C_ias","C_liq"), method="runsteady")$y
    
    # Integrate photosynthesis to leaf level
    leaf.vec[i] <- sum((pmin((k_c*X_c_z*out[,2])/(K_m+out[,2]),nx*(out[,2]*j_z/(4*out[,2]+8*Gamma)))-(pmin((k_c*X_c_z*out[,2])/(K_m+out[,2]),nx*(out[,2]*j_z/(4*out[,2]+8*Gamma)))*(Gamma/out[,2]))-Rd)*Vstrom*Sm)
    
    # Output intra-leaf CO2 consumption
    An.out[[i]] <- pmin((k_c*X_c_z*out[,2])/(K_m+out[,2]),nx*(out[,2]*j_z/(4*out[,2]+8*Gamma)))
    
    # Output intra-leaf CO2 concentration
    out.list[[i]] <- out
    
    # Output actual, potential, and maximum electron transport profiles
    j.out[[i]] <- data.frame(I_z,j_z,j_max_z) # write 
    
  }
  
  # Output chloroplast stroma volume
  out.Vstrom <- sum(Sm)*Vstrom 
  
  # Combine all output variables of interest into a list
  scenario.out[[j]] <- list(out.Vstrom,leaf.vec,out.list,An.out,j.out,Depth)
}

# Save list for each scenario into a list of lists
LeafOut.List[[n]] <- scenario.out
}

# Add leaf type + light type as names to list
names(LeafOut.List) = sapply(LIPNames, function(a) substring(a, 23,nchar(a)-4))

# Check which scenario minimizes the sum of squared errors with respect to observed data
sapply(1:length(LeafOut.List$Diffuse_Helianthus_Sun), function(i) sum((1e6*LeafOut.List$Diffuse_Helianthus_Sun[[i]][[2]]-Sun.Brodersen$A.dir)^2))
sapply(1:length(LeafOut.List$Direct_Helianthus_Shade), function(i) sum((1e6*LeafOut.List$Direct_Helianthus_Shade[[i]][[2]]-Shade.Brodersen$A.dif)^2))



###########################
# PLOT SIMULATION RESULTS #
###########################

# Plot leaf level light response curves
Cairo(file="Figs/Fig_LRC_wData_baseline.png", type='png', units='in', width=7.5, height=4, pointsize=12, dpi=300, bg='white')
par(mfrow=c(1,2),oma = c(3,4,0,0) + 0.1,mar = c(1,1,1,1) + 0.1)

# Sun leaf
plot(I.vec*1e6, LeafOut.List$Direct_Helianthus_Sun[[1]][[2]]*1e6, xlim=c(0,2100), ylim=c(-5,45), type = "l", xlab=c(expression(paste(Irradiance,' [','\u03BCmol ',m^2,' ',s^-1,']'))), ylab='', lwd=1.25, col='grey30')
points(d.mu[d.mu$trt=='sun_direct',]$I,d.mu[d.mu$trt=='sun_direct',]$A,col='grey30', bg='grey30', cex=1, pch=21)
lines(I.vec*1e6, LeafOut.List$Diffuse_Helianthus_Sun[[1]][[2]]*1e6, xlab=c(expression(paste(Irradiance,' [','\u03BCmol ',m^2,' ',s^-1,']'))), ylab='', col='grey70', lwd=1.25)
points(d.mu[d.mu$trt=='sun_diffuse',]$I,d.mu[d.mu$trt=='sun_diffuse',]$A,col='grey70', bg='grey70', cex=1, pch=21)
mtext(text=expression(paste(A[n],' [','\u03BCmol ',m^2,' ',s^-1,']')), line=1.5, side=2, cex=1.2, outer=TRUE)


# Shade leaf
plot(I.vec*1e6, LeafOut.List$Direct_Helianthus_Shade[[1]][[2]]*1e6, xlim=c(0,2100), ylim=c(-5,45), type = "l", xlab=c(expression(paste(Irradiance,' [','\u03BCmol ',m^2,' ',s^-1,']'))), ylab='', lwd=1.25, col='grey30', yaxt='n')
points(d.mu[d.mu$trt=='shade_direct',]$I,d.mu[d.mu$trt=='shade_direct',]$A, col='grey30', bg='grey30', cex=1, pch=21)
lines(I.vec*1e6, LeafOut.List$Diffuse_Helianthus_Shade[[1]][[2]]*1e6, xlab=c(expression(paste(Irradiance,' [','\u03BCmol ',m^2,' ',s^-1,']'))), ylab='', col='grey70', lwd=1.25, pch=21)
points(d.mu[d.mu$trt=='shade_diffuse',]$I,d.mu[d.mu$trt=='shade_diffuse',]$A,col='grey70', bg='grey70',cex=1, pch=21)
mtext(text=expression(paste(Illuminance,' [','\u03BCmol ',m^2,' ',s^-1,']')), line = 1.5, side=1, outer=TRUE)
dev.off()


# Plot CO2 liquid (C_liq) and intercellular airspace (C_ias) profiles
palette(rev(grey.colors(17, 0, 0.9)))
Cairo(file="Figs/FigS2_CO2_ProfilePanel_JmaxPropAbsorption.png", type='png', units='in', width=8, height=4.5, pointsize=15.5, dpi=150, bg='white')
par(mfrow=c(2,4),oma = c(5,4,0,0) + 0.1,mar = c(1,1,1,1) + 0.1)
plot(moltoppm*LeafOut.List$Direct_Helianthus_Sun[[1]][[3]][[1]][,2], LeafOut.List[[2]][[1]][[6]]*1e6, type = "l",xlim=c(100,350), ylim=c(0,275), col=1, ylab=c(expression(paste('distance from abaxial surface [','\u03BCm',']'))), xlab='', xaxt='n', lwd=1.25, axes=FALSE)
sapply(2:17, function(i) lines(moltoppm*LeafOut.List$Direct_Helianthus_Sun[[1]][[3]][[i]][,2], LeafOut.List[[2]][[1]][[6]]*1e6, type = "l", xlim=c(150,350), col=i, xlab='', ylab='', xaxt='n', yaxt='n', lwd=1.25))
axis(side=1, at = c(100,225,350,475), labels=FALSE)
axis(side=2, at=c(0,75,150,225,300))
box()
mtext(text=expression(paste('distance from abaxial surface [','\u03BCm',']')), line=-20,side=4)

plot(moltoppm*LeafOut.List$Diffuse_Helianthus_Sun[[1]][[3]][[1]][,2], LeafOut.List[[2]][[1]][[6]]*1e6, type = "l", xlim=c(100,350), ylim=c(0,275), col=1, ylab=c(expression(paste('distance from abaxial surface [','\u03BCm',']'))), xlab='', lwd=1.25, axes=FALSE)
sapply(2:17, function(i) lines(moltoppm*LeafOut.List$Diffuse_Helianthus_Sun[[1]][[3]][[i]][,2], LeafOut.List[[2]][[1]][[6]]*1e6, type = "l", xlim=c(150,350),col=i, xlab='', ylab='', axes=FALSE, lwd=1.25))
axis(side=1, at = c(100,225,350,475), labels=FALSE)
axis(side=2, at=c(0,75,150,225,300), labels=FALSE)
box()

plot(moltoppm*LeafOut.List$Direct_Helianthus_Shade[[1]][[3]][[1]][,2], LeafOut.List[[1]][[1]][[6]]*1e6, type = "l", xlim=c(100,350), ylim=c(0,275), col=1, ylab=c(expression(paste('distance from abaxial surface [','\u03BCm',']'))), xlab='', xaxt='n', yaxt='n',lwd=1.25, axes=FALSE)
sapply(2:17, function(i) lines(moltoppm*LeafOut.List$Direct_Helianthus_Shade[[1]][[3]][[i]][,2], LeafOut.List[[1]][[1]][[6]]*1e6, type = "l", xlim=c(150,350),col=i, xlab='', ylab='', xaxt='n', yaxt='n', lwd=1.25))
axis(side=1, at = c(100,225,350,475), labels=FALSE)
axis(side=2, at=c(0,75,150,225,300), labels=FALSE)
box()

plot(moltoppm*LeafOut.List$Diffuse_Helianthus_Shade[[1]][[3]][[1]][,2], LeafOut.List[[1]][[1]][[6]]*1e6, type = "l", xlim=c(100,350), ylim=c(0,275), col=1, ylab=c(expression(paste('distance from abaxial surface [','\u03BCm',']'))), xaxt='n', yaxt='n', xlab='', lwd=1.25, axes=FALSE)
sapply(2:17, function(i) lines(moltoppm*LeafOut.List$Diffuse_Helianthus_Shade[[1]][[3]][[i]][,2], LeafOut.List[[1]][[1]][[6]]*1e6, type = "l", xlim=c(150,350),col=i, xlab='', ylab='', xaxt='n', yaxt='n', lwd=1.25))
axis(side=1, at = c(100,225,350,475), labels=FALSE)
axis(side=2, at=c(0,75,150,225,300), labels=FALSE)
box()

# Plot CO2 intercellular airspace (IAS) profiles
plot(moltoppm*LeafOut.List$Direct_Helianthus_Sun[[1]][[3]][[1]][,1], LeafOut.List[[2]][[1]][[6]]*1e6, type = "l",xlim=c(100,350), ylim=c(0,275), col=1, ylab=c(expression(paste('distance from abaxial surface [','\u03BCm',']'))), xlab='', lwd=1.25, axes=FALSE)
sapply(2:17, function(i) lines(moltoppm*LeafOut.List$Direct_Helianthus_Sun[[1]][[3]][[i]][,1], LeafOut.List[[2]][[1]][[6]]*1e6, type = "l", xlim=c(150,350),col=i, xlab='', ylab='', xaxt='n', lwd=1.25))
axis(side=1, at = c(100,225,350,475))
axis(side=2, at=c(0,75,150,225,300))
box()

plot(moltoppm*LeafOut.List$Diffuse_Helianthus_Sun[[1]][[3]][[1]][,1], LeafOut.List[[2]][[1]][[6]]*1e6, type = "l", xlim=c(100,350), ylim=c(0,275), col=1, ylab=c(expression(paste('distance from abaxial surface [','\u03BCm',']'))), xlab='', yaxt='n', lwd=1.25, axes=FALSE)
sapply(2:17, function(i) lines(moltoppm*LeafOut.List$Diffuse_Helianthus_Sun[[1]][[3]][[i]][,1], LeafOut.List[[2]][[1]][[6]]*1e6, type = "l", xlim=c(150,350),col=i, xlab='', ylab='', yaxt='n', lwd=1.25))
axis(side=1, at = c(100,225,350,475))
axis(side=2, at=c(0,75,150,225,300), labels=FALSE)
box()

plot(moltoppm*LeafOut.List$Direct_Helianthus_Shade[[1]][[3]][[1]][,1], LeafOut.List[[1]][[1]][[6]]*1e6, type = "l", xlim=c(100,350), ylim=c(0,275), col=1, ylab=c(expression(paste('distance from abaxial surface [','\u03BCm',']'))), xlab='', yaxt='n',  lwd=1.25, axes=FALSE)
sapply(2:17, function(i) lines(moltoppm*LeafOut.List$Direct_Helianthus_Shade[[1]][[3]][[i]][,1], LeafOut.List[[1]][[1]][[6]]*1e6, type = "l", xlim=c(150,350),col=i, xlab='', ylab='', xaxt='n', lwd=1.25))
axis(side=1, at = c(100,225,350,475))
axis(side=2, at=c(0,75,150,225,300), labels=FALSE)
box()

plot(moltoppm*LeafOut.List$Diffuse_Helianthus_Shade[[1]][[3]][[1]][,1], LeafOut.List[[1]][[1]][[6]]*1e6, type = "l", xlim=c(100,350), ylim=c(0,275), col=1, ylab=c(expression(paste('distance from abaxial surface [','\u03BCm',']'))), yaxt='n', xlab='', lwd=1.25, axes=FALSE)
sapply(2:17, function(i) lines(moltoppm*LeafOut.List$Diffuse_Helianthus_Shade[[1]][[3]][[i]][,1], LeafOut.List[[1]][[1]][[6]]*1e6, type = "l", xlim=c(150,350),col=i, xlab='', ylab='', xaxt='n', yaxt='n', lwd=1.25))
axis(side=1, at = c(100,225,350,475))
axis(side=2, at=c(0,75,150,225,300), labels=FALSE)
box()
mtext(text=expression(paste(CO[2],' concentration [ppm]')), line=2, side=1, outer=TRUE, cex=0.8)
mtext(text=expression(paste('distance from abaxial surface [','\u03BCm',']')), line=2, side=2, outer=TRUE, cex=0.8)
dev.off()

# Plot the difference between direct and diffuse light scenarios for the sun and shade leaves
Cairo(file="Figs/TEST.png", type='png', units='in', width=7.5, height=4, pointsize=12, dpi=300, bg='white')
par(mfrow=c(1,2), oma = c(3,3,0,0) + 0.1, mar = c(1,1,1,1) + 0.1)
plot(moltoppm*(LeafOut.List$Direct_Helianthus_Sun[[1]][[3]][[1]][,2]-LeafOut.List$Diffuse_Helianthus_Sun[[1]][[3]][[1]][,2]), LeafOut.List[[2]][[1]][[6]]*1e6, type = "l",xlim=c(-60,30), ylim=c(0,275), col=1, ylab=c(expression(paste('distance from abaxial surface [','\u03BCm',']'))), xlab=c(expression(paste(dC[liq], ' [ppm]'))), lwd=1.5)
sapply(2:17, function(i) lines(moltoppm*(LeafOut.List$Direct_Helianthus_Sun[[1]][[3]][[i]][,2]-LeafOut.List$Diffuse_Helianthus_Sun[[1]][[3]][[i]][,2]), LeafOut.List[[2]][[1]][[6]]*1e6, type = "l", xlim=c(-60,30),col=i, xlab='', ylab='', xaxt='n', yaxt='n', lwd=1.5))
text(x=-40,y=250,labels='sun leaf', cex=1)

plot(moltoppm*(LeafOut.List$Direct_Helianthus_Shade[[1]][[3]][[1]][,2]-LeafOut.List$Diffuse_Helianthus_Shade[[1]][[3]][[1]][,2]), LeafOut.List[[1]][[1]][[6]]*1e6, type = "l",xlim=c(-60,30), ylim=c(0,275), col=1, ylab='', xlab=c(expression(paste(dC[liq], ' [ppm]'))), lwd=1.5, yaxt='n')
sapply(2:17, function(i) lines(moltoppm*(LeafOut.List$Direct_Helianthus_Shade[[1]][[3]][[i]][,2]-LeafOut.List$Diffuse_Helianthus_Shade[[1]][[3]][[i]][,2]), LeafOut.List[[1]][[1]][[6]]*1e6, type = "l", xlim=c(-60,30),col=i, xlab='', ylab='', xaxt='n', yaxt='n', lwd=1.5))
text(x=-40,y=250,labels='shade leaf', cex=1)
mtext(c(expression(paste('distance from abaxial surface [','\u03BCm',']'))), side = 2, outer = TRUE, cex = 1, line = 1.5,col = "grey20")
mtext(expression(paste(dC[liq], ' [ppm]')), side = 1, outer = TRUE, cex = 1, line = 1.5,col = "grey20")
dev.off()

# Plot parametric sensitivity analysis
# Estimate theta and Amax
#scenarios.df$A750 = 0
An_avg = 0
scenarios.df$Theta = 0
scenarios.df$phiCO2 = 0
A_I <- function(I,theta,phi_co2) (phi_co2*I + Amax - sqrt((phi_co2*I+Amax)^2-4*theta*phi_co2*I*Amax))/(2*theta)
for(j in 1:length(LeafOut.List)){
  for(i in 1:nrow(scenarios.df)){
    An = LeafOut.List[[j]][[i]][[2]]*1e6
    I = I.vec*1e6
    Amax = max(LeafOut.List[[j]][[i]][[2]])*1e6
    #LeafOut.List[[j]]$A750[i] = An[I(which(I==750))]
    LeafOut.List[[j]]$An_avg[i] = trapz(I[I <= 1500], An[I <= 1500]) / max(I[I <= 1500])
    LeafOut.List[[j]]$Amax[i] = Amax
    df = data.frame(I,An)
    LeafOut.List[[j]]$Theta[i] <- coef(nls(An ~ A_I(I, theta, phi_co2), data=df, start=list(theta=0.999, phi_co2=0.04), lower=list(theta=0), upper=list(theta=1), trace=T, alg='port'))[[1]]
    LeafOut.List[[j]]$phiCO2[i] <- coef(nls(An ~ A_I(I, theta, phi_co2), data=df, start=list(theta=0.999, phi_co2=0.04), lower=list(theta=0), upper=list(theta=1), trace=T, alg='port'))[[2]]
  }
}


# Plot parametric sensitivity analysis figure
scen.sun <- data.frame(scenarios.df$PARM[1:31],(LeafOut.List[[2]]$An_avg-LeafOut.List[[4]]$An_avg))[-c(30,31),]
scen.sun$Theta <- (LeafOut.List[[2]]$Theta-LeafOut.List[[4]]$Theta)[-c(30,31)]
scen.shade <- data.frame(scenarios.df$PARM[1:31],(LeafOut.List[[1]]$An_avg-LeafOut.List[[3]]$An_avg))[-c(30,31),]
scen.shade$Theta <- (LeafOut.List[[1]]$Theta-LeafOut.List[[3]]$Theta)[-c(30,31)]
colnames(scen.sun) = colnames(scen.shade) = c('PARM','dAn_avg','Theta')
scen.sun$HL = scen.shade$HL = c('baseline',rep(c('low','high'),14))
var.reorder = c('.Alpha','.Beta','.Jmax','.k_c','.K_m','.X_c','.frac_pal','.g_liq','.por_pal','.por_spg','.Sm_pal','.Sm_spg','.tort','.Vstrom')

p1 <- ggplot(scen.sun[-1,], aes(x=PARM,y=dAn_avg)) +  geom_hline(aes(yintercept=scen.sun[1,]$dAn_avg), linetype='dotted') + geom_line(size=3, col='grey80') + geom_point(aes(x=PARM,y=dAn_avg, col=HL), size=3)  + xlab('sun leaf') + theme_bw() + scale_x_discrete(labels=c('\u03B1','\u03B2',bquote(J[max]),bquote(k[c]),bquote(K[m]),bquote(X[c]),bquote(f[pal]),bquote(g[liq]),bquote('\u03C6'[pal]),bquote('\u03C6'[spg]),bquote(S[m][','][pal]),bquote(S[m][','][spg]),bquote('\u03c4'),bquote(V[strom])), limits=var.reorder) + ylim(-4,1)  + scale_color_manual(values=c("black","grey60")) + theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=13)) + ylab('') + coord_flip()

p2 <- ggplot(scen.shade[-1,], aes(x=PARM,y=dAn_avg)) +  geom_hline(aes(yintercept=scen.shade[1,]$dAn_avg), linetype='dotted') + geom_line(size=3, col='grey80') + geom_point(aes(x=PARM,y=dAn_avg, col=HL), size=3) + ylab(expression(paste(bar(A)['n,diffuse'],' - ',bar(A)['n,direct'],' [','\u03BCmol ',m^2,' ',s^-1,']'))) + xlab('shade leaf') + theme_bw() + scale_x_discrete(labels=c('\u03B1','\u03B2',bquote(J[max]),bquote(k[c]),bquote(K[m]),bquote(X[c]),bquote(f[pal]),bquote(g[liq]),bquote('\u03C6'[pal]),bquote('\u03C6'[spg]),bquote(S[m][','][pal]),bquote(S[m][','][spg]),bquote('\u03c4'),bquote(V[strom])), limits=var.reorder) + coord_flip() + ylim(-4,1) + scale_color_manual(values=c("black","grey60")) + theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=13))

p3 <- ggplot(scen.sun[-1,], aes(x=PARM,y=Theta)) +  geom_hline(aes(yintercept=scen.sun[1,]$Theta), linetype='dotted') + geom_line(size=3, col='grey80') + geom_point(aes(x=PARM,y=Theta, col=HL), size=3)  + xlab('') + theme_bw() + scale_x_discrete(labels=c('\u03B1','\u03B2',bquote(J[max]),bquote(k[c]),bquote(K[m]),bquote(X[c]),bquote(f[pal]),bquote(g[liq]),bquote('\u03C6'[pal]),bquote('\u03C6'[spg]),bquote(S[m][','][pal]),bquote(S[m][','][spg]),bquote('\u03c4'),bquote(V[strom])), limits=var.reorder) + coord_flip() + ylim(-0.1,0)  + scale_color_manual(values=c("black","grey60")) + theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=13)) + ylab('') + annotate("text", y=-3.25, x=13, label="sun leaf", size=4.5)

p4 <- ggplot(scen.shade[-1,], aes(x=PARM,y=Theta)) +  geom_hline(aes(yintercept=scen.shade[1,]$Theta), linetype='dotted') + geom_line(size=3, col='grey80') + geom_point(aes(x=PARM,y=Theta, col=HL), size=3) + ylab(expression(paste('\u03f4'[diffuse],' - ','\u03f4'[direct]))) + xlab('') + theme_bw() + scale_x_discrete(labels=c('\u03B1','\u03B2',bquote(J[max]),bquote(k[c]),bquote(K[m]),bquote(X[c]),bquote(f[pal]),bquote(g[liq]),bquote('\u03C6'[pal]),bquote('\u03C6'[spg]),bquote(S[m][','][pal]),bquote(S[m][','][spg]),bquote('\u03c4'),bquote(V[strom])), limits=var.reorder) + coord_flip() + ylim(-0.1,0) + scale_color_manual(values=c("black","grey60")) + theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size=13)) + annotate("text", y=-3, x=13, label="shade leaf", size=4.5)

grid.arrange(p1,p3,p2,p4,ncol=2)

# Save parameter sensitivity analysis figure
png(file="../../Figs/SensitivityAnalysis_updated.png", units='px', width=1100, height=1000, pointsize=12, res=150, bg='white')
grid.arrange(p1,p3,p2,p4,ncol=2)
dev.off()

# Examine interaction between maximum e- transport rate (j_max), realized e- rate (j_e), and fraction of light absorbed by PSII that drives e- transport (I_z)
png(file="Figs/FigS7a_ElectronTransport_ConstantPhiPSII_SunLeaf.png", units='px', width=1100, height=600, pointsize=12, res=150, bg='white')
png(file="Figs/Fig1_ConceptualDiagram.png", units='px', width=700, height=750, pointsize=12, res=150, bg='white')

I.sub = c(5,13,15,17) #select subset of irradiance values
Depth = LeafOut.List[[2]][[1]][[6]]
par(mfrow=c(2,4),oma = c(4,4,2,2) + 0.1,mar = c(1,1,1,1) + 0.1,cex.lab=1.5, cex.axis=1.25)
for(j in c('Direct_Helianthus_Sun','Diffuse_Helianthus_Sun')) {
  for(i in 1:4){
    j.df = sapply(I.sub, function(i) unlist(LeafOut.List[[which(names(LeafOut.List)==j)]][[1]][[5]][[i]]))
    Iz = nx*j.df[1:nx,]; je = nx*j.df[(nx+1):(2*nx),]; jmax = nx*j.df[(2*nx+1):(3*nx),]
    if(i == 1 & j != 'Diffuse_Helianthus_Sun') {
      plot(jmax[,i], Depth*1e6, xlim=c(0,15), ylim=c(0,275), type='l', lwd=3., xlab='', ylab='', xaxt='n')
      lines(Iz[,i], Depth*1e6, col='gold', lwd=3.)
      lines(je[,i],Depth*1e6, col='red', lwd=3.)
    }
    if(i != 1 & j != 'Diffuse_Helianthus_Sun') {
      plot(jmax[,i],Depth*1e6, xlim=c(0,15),  ylim=c(0,275), type='l', lwd=3., xlab='', ylab='', xaxt='n', yaxt='n')
      lines(Iz[,i],Depth*1e6, col='gold', lwd=3.)
      lines(je[,i],Depth*1e6, col='red', lwd=3.) 
    }
    if(i == 1 & j == 'Diffuse_Helianthus_Sun') {
      plot(jmax[,i],Depth*1e6, xlim=c(0,15), ylim=c(0,275), type='l', lwd=3., xlab='', ylab='')
      lines(Iz[,i],Depth*1e6, col='gold', lwd=3.)
      lines(je[,i],Depth*1e6, col='red', lwd=3.)
      legend(5.5,130,legend = c(expression(paste(I[e])),expression(paste(j[e])),expression(paste(j[max]))), col=c('gold','red','black'),lty=1, bty='n', cex=1.35, lwd=1.5, y.intersp=.9)
    }
    if(i != 1 & j == 'Diffuse_Helianthus_Sun') {
      plot(jmax[,i],Depth*1e6, xlim=c(0,15), ylim=c(0,275), type='l', lwd=3., xlab='', ylab='', yaxt='n')
      lines(Iz[,i],Depth*1e6, col='gold', lwd=3.)
      lines(je[,i],Depth*1e6, col='red', lwd=3.) 
    }
  }
}
mtext(c(expression(paste('distance from abaxial surface [','\u03BCm',']'))), side = 2, outer = TRUE, cex = 1.2, line = 1.5,col = "grey20")
mtext(c(expression(paste('mol e- or photons ',m^-3,' ',s^-1))), side = 1, outer = TRUE, cex = 1.2, line = 2,col = "grey20")
mtext('2000 PPFD', side = 3, outer = TRUE, cex = 0.8, line = 0, at = c(0.875,0.875), col = "grey20")
mtext('1250 PPFD', side = 3, outer = TRUE, cex = 0.8, line = 0, at = c(0.625,0.625), col = "grey20")
mtext('750 PPFD', side = 3, outer = TRUE, cex = 0.8, line = 0, at = c(0.375,0.375), col = "grey20")
mtext('250 PPFD', side = 3, outer = TRUE, cex = 0.8, line = 0, at = c(0.125,0.125), col = "grey20")
dev.off()

png(file="Figs/FigS7b_ElectronTransport_ConstantPhiPSII_ShadeLeaf.png", units='px', width=1100, height=600, pointsize=12, res=150, bg='white')
I.sub = c(5,13,15,17) #select subset of irradiance values
Depth = LeafOut.List[[1]][[1]][[6]]
par(mfrow=c(2,4),oma = c(4,4,2,2) + 0.1,mar = c(1,1,1,1) + 0.1,cex.lab=1.5, cex.axis=1.25)
for(j in c('Direct_Helianthus_Shade','Diffuse_Helianthus_Shade')) {
  for(i in 1:4){
    j.df = sapply(I.sub, function(i) unlist(LeafOut.List[[which(names(LeafOut.List)==j)]][[1]][[5]][[i]]))
    Iz = j.df[1:nx,]; je = j.df[(nx+1):(2*nx),]; jmax = j.df[(2*nx+1):(3*nx),]
    if(i == 1 & j != 'Diffuse_Helianthus_Shade') {
      plot(jmax[,i], Depth*1e6, xlim=c(0,0.06), ylim=c(0,275), type='l', lwd=3., xlab='', ylab='', xaxt='n')
      lines(Iz[,i], Depth*1e6, col='gold', lwd=3.)
      lines(je[,i],Depth*1e6, col='red', lwd=3.)
    }
    if(i != 1 & j != 'Diffuse_Helianthus_Shade') {
      plot(jmax[,i],Depth*1e6, xlim=c(0,0.06), ylim=c(0,275), type='l', lwd=3., xlab='', ylab='', xaxt='n', yaxt='n')
      lines(Iz[,i],Depth*1e6, col='gold', lwd=3.)
      lines(je[,i],Depth*1e6, col='red', lwd=3.) 
    }
    if(i == 1 & j == 'Diffuse_Helianthus_Shade') {
      plot(jmax[,i],Depth*1e6, xlim=c(0,0.06), ylim=c(0,275), type='l', lwd=3., xlab='', ylab='')
      lines(Iz[,i],Depth*1e6, col='gold', lwd=3.)
      lines(je[,i],Depth*1e6, col='red', lwd=3.)
      legend(0.023,140,legend = c(expression(paste(I[e])),expression(paste(j[e])),expression(paste(j[max]))), col=c('gold','red','black'),lty=1, bty='n', cex=1.35, lwd=1.5, y.intersp=.9)
    }
    if(i != 1 & j == 'Diffuse_Helianthus_Shade') {
      plot(jmax[,i],Depth*1e6, xlim=c(0,0.06), ylim=c(0,275), type='l', lwd=3., xlab='', ylab='', yaxt='n')
      lines(Iz[,i],Depth*1e6, col='gold', lwd=3.)
      lines(je[,i],Depth*1e6, col='red', lwd=3.) 
    }
  }
}
mtext(c(expression(paste('distance from abaxial surface [','\u03BCm',']'))), side = 2, outer = TRUE, cex = 1, line = 1.5,col = "grey20")
mtext(c(expression(paste('\u03BCmol e- or photons ',m^-3,' ',s^-1))), side = 1, outer = TRUE, cex = 1, line = 2,col = "grey20")
mtext('2000 PPFD', side = 3, outer = TRUE, cex = 0.8, line = 0, at = c(0.875,0.875), col = "grey20")
mtext('1250 PPFD', side = 3, outer = TRUE, cex = 0.8, line = 0, at = c(0.625,0.625), col = "grey20")
mtext('750 PPFD', side = 3, outer = TRUE, cex = 0.8, line = 0, at = c(0.375,0.375), col = "grey20")
mtext('250 PPFD', side = 3, outer = TRUE, cex = 0.8, line = 0, at = c(0.125,0.125), col = "grey20")
dev.off()

# Plot Jmax scenarios
# Jmax proportional to absorption profile

# Sun leaf
png(file="Figs/Fig11_LRC_Cliq_JmaxPropChlDistr.png", units='px', width=1000, height=550, pointsize=12, res=150, bg='white')
par(mfrow=c(1,2),oma = c(4,4,2,4) + 0.1,mar = c(1,1,1,1) + 0.1)

#layout(matrix(c(1,2,1,3), 2, 2, byrow=TRUE))
#par(oma = c(2,2,2,2))
plot(I.vec*1e6, LeafOut.List$Direct_Helianthus_Sun[[1]][[2]]*1e6, xlim=c(0,2100), ylim=c(-5,45), type = "l", xlab='', ylab='', lwd=0.75, col='grey30', lty='dashed', xaxt='n')
axis(side=1, at=c(0,1000,2000))
points(d.mu[d.mu$trt=='sun_direct',]$I,d.mu[d.mu$trt=='sun_direct',]$A,col='grey30', bg='grey30', cex=1, pch=21)
lines(I.vec*1e6, LeafOut.List$Diffuse_Helianthus_Sun[[1]][[2]]*1e6, xlab='', ylab='', col='grey70', lwd=0.75, lty='dashed', xaxt='n')
points(d.mu[d.mu$trt=='sun_diffuse',]$I,d.mu[d.mu$trt=='sun_diffuse',]$A,col='grey70', bg='grey70', cex=1, pch=21)
mtext(text=expression(paste(A[n],' [','\u03BCmol ',m^2,' ',s^-1,']')), line=1.5, side=2, cex=1, outer=TRUE)
#mtext(text=expression(paste(Irradiance,' [','\u03BCmol ',m^2,' ',s^-1,']')), line=1.5, side=1, cex=1.2, outer=TRUE)
mtext(text=expression(paste(Irradiance,' [','\u03BCmol ',m^2,' ',s^-1,']')), line=1.5, side=1, cex=1, outer=TRUE, adj=0.12)

plot(moltoppm*(LeafOut.List$Direct_Helianthus_Sun[[1]][[3]][[1]][,2]-LeafOut.List$Diffuse_Helianthus_Sun[[1]][[3]][[1]][,2]), Depth*1e6, type = "l",xlim=c(-60,30), col=1, ylab=c(expression(paste('distance from abaxial surface [','\u03BCm',']'))), xlab=c(expression(paste(dC[liq], ' [ppm]'))), lwd=1.5, yaxt='n')
sapply(2:17, function(i) lines(moltoppm*(LeafOut.List$Direct_Helianthus_Sun[[1]][[3]][[i]][,2]-LeafOut.List$Diffuse_Helianthus_Sun[[1]][[3]][[i]][,2]), Depth*1e6, type = "l", xlim=c(-60,30),col=i, xlab='', ylab='', xaxt='n', yaxt='n', lwd=1.5))
axis(side=4, at = c(0,75,150,225))
mtext(text=expression(paste(dC[liq],' [ppm]')), line=1.5, side=1, cex=1, outer=TRUE, adj=0.8)
mtext(text=expression(paste('distance from abaxial surface [','\u03BCm',']')), line=1.5, side=4, cex=1, outer=TRUE)
dev.off()

#############################
# Canopy-level Calculations #
#############################

# Define extinction coefficients
K_dir = 1.06*0.25 # Direct light extinction coefficient
K_dif = 0.82*0.25 # Diffuse light extinction coefficient

# Define sequence of cumulative LAI values for simulation
L = seq(0, 5, 0.1)

# Define function for calculating irradiance at different cumulative LAIs within canopy
I_L = function(L, I_0, K) I_0 * exp(-K * L)

# Calculate irradiance at different cumulative LAIs within canopy
I_L_canopy_dir = sapply(L, I_L, I_0 = 1500, K = K_dir)
I_L_canopy_dif = sapply(L, I_L, I_0 = 1500, K = K_dif)

# Define assimilation function
A_I = function(I, theta, phi_co2, Amax) (phi_co2*I + Amax - sqrt((phi_co2*I+Amax)^2-4*theta*phi_co2*I*Amax))/(2*theta)

# Direct sun scenario
leaf_parms = c(LeafOut.List[[4]]$Theta[1], LeafOut.List[[4]]$phiCO2[1], LeafOut.List[[4]]$Amax[1])
A_can_DirSun = sapply(I_L_canopy_dir, A_I, theta = leaf_parms[1], phi_co2 = leaf_parms[2], Amax = leaf_parms[3])

# Diffuse sun scenario
leaf_parms = c(LeafOut.List[[2]]$Theta[1], LeafOut.List[[2]]$phiCO2[1], LeafOut.List[[2]]$Amax[1])
A_can_DifSun = sapply(I_L_canopy_dif, A_I, theta = leaf_parms[1], phi_co2 = leaf_parms[2], Amax = leaf_parms[3])

# Direct shade scenario
leaf_parms = c(LeafOut.List[[3]]$Theta[1], LeafOut.List[[3]]$phiCO2[1], LeafOut.List[[3]]$Amax[1])
A_can_DirShd = sapply(I_L_canopy_dir, A_I, theta = leaf_parms[1], phi_co2 = leaf_parms[2], Amax = leaf_parms[3])

# Diffuse shade scenario
leaf_parms = c(LeafOut.List[[1]]$Theta[1], LeafOut.List[[1]]$phiCO2[1], LeafOut.List[[3]]$Amax[1])
A_can_DifShd = sapply(I_L_canopy_dif, A_I, theta = leaf_parms[1], phi_co2 = leaf_parms[2], Amax = leaf_parms[3])

# Calculate cross-over point between diffuse and direct light
Sun_cross = which(abs(A_can_DirSun - A_can_DifSun) == min(abs(A_can_DirSun - A_can_DifSun)))
Shd_cross = which(abs(A_can_DirShd - A_can_DifShd) == min(abs(A_can_DirShd - A_can_DifShd)))

##
# Plot output
##

# Plot canopy light availability curves
png('CanopyLevel_LowExtinction.png', width = 7.5, height = 2.25, units = 'in', pointsize = 7, res = 150)
par(mfrow = c(1,3), oma = c(1, 1, 0, 0), mar = c(3, 4, 2, 2), cex = 1.15)
plot(L, I_L_canopy_dir, type = 'l', xlab = '', ylab = '', xlim = c(0, 5), ylim = c(0, 1500))
lines(L, I_L_canopy_dif, lty = 2)
mtext('', side = 1, line = 3)
mtext(expression(paste('I [', mu, 'mol ', ' ', m^-2, s^-1, ']')), side = 2, line = 2.5, cex = 1.2)
legend(x = 2.5,
       y = 1600,
       lty = c(1, 2),
       legend = c('direct', 'diffuse'),
       bty = 'n',
       x.intersp = 0.3,
       y.intersp = 1.5,
       xjust = 0,
       seg.len = 2.5
)

# Plot canopy light availability curves
plot(L[1:25], A_can_DirSun[1:25], type = 'l', xlab = '', ylab = '', ylim = c(0,35))
lines(L[1:25], A_can_DifSun[1:25], lty = 2)
polygon(x = c(L[Sun_cross], L[Sun_cross], max(L), max(L)),
        y = c(0, 35, 35, 0),
        col = rgb(0, 0, 0,0.2),
        border = NA)
mtext(expression(paste(A[n], ' [', mu, 'mol ', m^-2, ' ', s^-1, ']')), side = 2, line = 2.5, cex = 1.2)
text(x = 3, y = 33, labels = 'sun-grown leaves only')

# Plot canopy light availability curves
plot(L, A_can_DirShd, type = 'l', xlab = '', ylab = '', ylim = c(0,35))
lines(L, A_can_DifShd, lty = 2)
mtext(expression(paste(A[n], ' [', mu, 'mol ' , ' ', m^-2, s^-1, ']')), side = 2, line = 2.5, cex = 1.2)
text(x = 3, y = 33, labels = 'shade-grown leaves only')
mtext(expression(paste('cumulative LAI [', m^2, ' ', m^-2, ']')), side = 1, line = -0.5, at = 0.52, outer = TRUE, cex = 1.2)

dev.off()