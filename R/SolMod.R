#################################################################################
# These functions are written by D. Yang
# Singapore Institute of Manufacturing Technology (SIMTech)
# emails: yangdazhi.nus@gmail.com; yangdz@SIMTech.a-star.edu.sg
#################################################################################

#################################################################################
# function to calculate sun position, extraterrestial irradiance,
# clear sky (Perez-Ineichen) irradiance and apparent solar time
#################################################################################
solpos <- function(Tm, lat, lon, tz, tilt, orientation, LT, alt = 0)
{
	#tilt angle between 0 and 90
	#orientation between 0 and 360, north 0, east positive(clockwise)
	jd <- JD(Tm)
	sunv <- sunvector(jd, lat, lon, tz)
	azi <- round(sunpos(sunv)[,1],3) #azimuth of the sun
	zen <- round(sunpos(sunv)[,2],3) #zenith angle
	surface.norm <- normalvector(tilt, orientation)
	inc <- round(as.numeric(degrees(acos(sunv%*% as.vector(surface.norm)))),3)
	dec <- declination(jd)*pi/180
	re <- 1.000110+0.034221*cos(dec)+0.001280*sin(dec)+0.000719*cos(2*dec)+0.000077*sin(2*dec)
	Io <- round(1362*re, 3) #extraterrestrial direct normal irradiance
	Ioh <- round(1362*re*cos(zen*pi/180), 3) #horizontal extraterrestrial irradiance
	Ioh <- ifelse(zen>=90, 0, Ioh)

  # Equation of time (L. O. Lamm, 1981, Solar Energy 26, p465)
  dn <- round(as.numeric(format(Tm, "%j"))-1 + as.numeric(format(Tm, "%H"))/24, 3)
  coef <- matrix(c(0, 0.00020870, 0,
            			 1, 0.0092869, -0.12229,
  	    			     2, -0.052258, -0.15698,
  			      		 3, -0.0013077, -0.0051602,
  					       4, -0.0021867, -0.0029823,
  					       5, -0.00015100, -0.00023463), ncol = 3, byrow = TRUE)
  EOT <- rowSums(sapply(1:6, function(i) coef[i,2]*cos(2*pi*coef[i,1]*dn/365.25) + coef[i,3]*sin(2*pi*coef[i,1]*dn/365.25)))*60 #EOT in minutes
  ast <- Tm - 4*60*(tz*15-lon) + EOT*60 #apparent solar time in POSIXct format

	#Perez-Ineichen clear sky model (Ineichen and Perez, 2002), monthly Linke turbidity can be obtained from SoDa service (following: Gueymard and Ruiz-Aries, 2016)
	fh1 <- exp(-alt/8000)
	fh2 <- exp(-alt/1250)
	cg1 <- (0.0000509*alt + 0.868)
	cg2 <- (0.0000392*alt + 0.0387)
	AM <- 1/(cos(zen*pi/180)+0.50572*(96.07995 - zen)^(-1.6364))
	Ics <- cg1*Io*cos(zen*pi/180)*exp(-cg2*AM*(fh1 + fh2*(LT-1)))*exp(0.01*AM^1.8)
	Ics <- ifelse(zen>=90, 0, round(Ics, 3))
	Icsd <- (0.664+0.163/fh1)*Io*exp(-0.09*(LT-1)*AM)*cos(zen*pi/180)
  Icsd <- ifelse(zen>=90, 0, round(Icsd, 3))

	out = list(zen, inc, Io, Ioh, azi, Ics, Icsd, ast)
	names(out) = c("zenith", "incidence", "Io", "Ioh", "azimuth", "Ics", "Icsd", "ast")
	return(out)
}

#degree to radian
d2r <-function(x)
{
  x*pi/180
}

#radian to degree
r2d <-function(x)
{
  x*180/pi
}

#################################################################################
# function to load the original set of parameters fitted by Perez et.al. (1990)
#################################################################################
load.Perez.par <- function(model = "Perez90")
{
  if (model == "Perez90") #Simplified model using data from 10 American and 3 European sites
  {
    F.mat = matrix(c(
      -0.008,0.588,-0.062,-0.060,0.072,-0.022,
      0.130,0.683,-0.151,-0.019,0.066,-0.029,
      0.330,0.487,-0.221,0.055,-0.064,-0.026,
      0.568,0.187,-0.295,0.109,-0.152,-0.014,
      0.873,-0.392,-0.362,0.226,-0.462,0.001,
      1.133,-1.237,-0.412,0.288,-0.823,0.056,
      1.060,-1.600,-0.359,0.264,-1.127,0.131,
      0.678,-0.327,-0.250,0.156,-1.377,0.251
    ), ncol = 6, byrow = TRUE)

    epsilon.mat = matrix(c(1, 1.065, 1.065, 1.23, 1.23, 1.5, 1.5, 1.95, 1.95, 2.8, 2.8, 4.5, 4.5, 6.2, 6.2, 12), ncol = 2, byrow = TRUE)
  }

  if (model == "Perez87a") #Simplified circumsolar model with 25 half-angle circumsolar region, using 2 French sites
  {
    F.mat = matrix(c(
      -0.011, 0.748, -0.080, -0.048, 0.073, -0.024,
      -0.038, 1.115, -0.109, -0.023, 0.106, -0.037,
      0.166, 0.909, -0.179, 0.062, -0.021, -0.050,
      0.419, 0.646, -0.262, 0.140, -0.167, -0.042,
      0.710, 0.025, -0.290, 0.243, -0.511, -0.004,
      0.857, -0.370, -0.279, 0.267, -0.792, 0.076,
      0.743, -0.073, -0.228, 0.231, -1.180, 0.199,
      0.421, -0.661, 0.097, 0.119, -2.125, 0.446
    ), ncol = 6, byrow = TRUE)

    epsilon.mat = matrix(c(1, 1.056, 1.056, 1.253, 1.253, 1.586, 1.586, 2.134, 2.134, 3.230, 3.230, 5.980, 5.980, 10.080, 10.080, 10000), ncol = 2, byrow = TRUE)
  }

  if (model == "Perez87b") #Simplified model using 2 French sites
  {
    F.mat = matrix(c(
      0.041, 0.621, -0.105, -0.040, 0.074, -0.031,
      0.054, 0.966, -0.166, -0.016, 0.114, -0.045,
      0.227, 0.866, -0.250, 0.069, -0.002, -0.062,
      0.486, 0.670, -0.373, 0.148, -0.137, -0.056,
      0.819, 0.106, -0.465, 0.268, -0.497, -0.029,
      1.020, -0.260, -0.514, 0.306, -0.804, 0.046,
      1.009, -0.708, -0.433, 0.287, -1.286, 0.166,
      0.936, -1.121, -0.352, 0.226, -2.449, 0.383
    ), ncol = 6, byrow = TRUE)

    epsilon.mat = matrix(c(1, 1.056, 1.056, 1.253, 1.253, 1.586, 1.586, 2.134, 2.134, 3.230, 3.230, 5.980, 5.980, 10.080, 10.080, 10000), ncol = 2, byrow = TRUE)
  }

  list(F.mat = F.mat, epsilon.mat = epsilon.mat)
}


#################################################################################
# Liu and Jordan 1961 (s)
#################################################################################
Liu <- function(data)
{
  s <- d2r(data$s)
  Rd <- (1+cos(s))*0.5
  Rd
}

#################################################################################
# Bugler 1977 (s, theta, z)
#################################################################################
Bugler1 <- function(data)
{
  theta <- d2r(data$theta)
  z <- d2r(data$Z)
  s <- d2r(data$s)
  Ih <- data$Gh - data$Dh
  rb <- pmax(0, cos(theta)/cos(z))
  Rd <- (cos(s/2))^2 + 0.05*cos(theta)*Ih*rb/data$Dh
  Rd
}

Bugler2 <- function(data)
{
  theta <- d2r(data$theta)
  z <- d2r(data$Z)
  s <- d2r(data$s)
  Ih <- data$Gh - data$Dh
  rb <- pmax(0, cos(theta)/cos(z))
  Rd <- (1-0.05*Ih/data$Dh)*(cos(s/2))^2 + 0.05*cos(theta)*Ih*rb/data$Dh
  Rd
}

#################################################################################
# Temps and Coulson 1977 (s, theta, z)
#################################################################################
Temps <- function(data)
{
	theta <- d2r(data$theta)
	z <- d2r(data$Z)
	s <- d2r(data$s)
	Rd <- (cos(s/2))^2*(1+(sin(s/2))^3)*(1+(cos(theta))^2*(sin(z))^3)
	Rd
}

#################################################################################
# Klucher 1979 (s, theta, z, Gh, Dh)
#################################################################################
Klucher <- function(data)
{
	Ff <- 1-(data$Dh/data$Gh)^2
	theta <- d2r(data$theta)
	z <- d2r(data$Z)
	s <- d2r(data$s)
	Rd <- (cos(s/2))^2*(1+Ff*(sin(s/2))^3)*(1+Ff*(cos(theta))^2*(sin(z))^3)
	Rd
}

#################################################################################
# Hay and Davies 1980, 1993 (s, theta, z, Ih, Ioh)
#################################################################################
Hay1 <- function(data)
{
	AI <- (data$Gh - data$Dh)/data$Ioh
	theta <- d2r(data$theta)
	z <- d2r(data$Z)
	s <- d2r(data$s)
	rb <- pmax(0, cos(theta)/cos(z))
	Rd <- (1-AI)*(cos(s/2))^2 + AI*rb
	Rd
}

Hay2 <- function(data)
{
	AI <- data$I/1362
	theta <- d2r(data$theta)
	z <- d2r(data$Z)
	s <- d2r(data$s)
	rb <- pmax(0, cos(theta)/cos(z))
	Rd <- (1-AI)*(cos(s/2))^2 + AI*rb
	Rd
}

#################################################################################
# Steven and Unsworth 1979, 1980 (s, theta, z)
#################################################################################
Steven1 <- function(data)
{
	theta <- d2r(data$theta)
	z <- d2r(data$Z)
	s <- d2r(data$s)
	rb <- pmax(0, cos(theta)/cos(z))
	steven <- sin(s) - s*cos(s) -pi*(sin(s/2))^2
	Rd <- 0.51*rb + (1-0.51)*((cos(s/2))^2 - 0.4395708*steven)
	Rd
}

Steven2 <- function(data)
{
	theta <- d2r(data$theta)
	z <- d2r(data$Z)
	s <- d2r(data$s)
	rb <- pmax(0, cos(theta)/cos(z))
	steven <- sin(s) - s*cos(s) -pi*(sin(s/2))^2
	#select a and b values
	ab.range <- d2r(c(35,45,55,65))
	tmp <- abs(z - matrix(rep(ab.range, length(z)), ncol=4, byrow = TRUE))
	choice <- apply(tmp, 1, which.min)
	a <- c(0.63, 0.6, 0.53, 0.46); b <- c(-1.04, -1, -0.9, -0.85);
	Rd <- a[choice]*rb + (1-a[choice])*((cos(s/2))^2 + 2*b[choice]/pi/(3 + 2*b[choice])*steven)
	Rd
}

Steven3 <- function(data)
{
	theta <- d2r(data$theta)
	z <- d2r(data$Z)
	s <- d2r(data$s)
	rb <- pmax(0, cos(theta)/cos(z))
	steven <- sin(s) - s*cos(s) -pi*(sin(s/2))^2
	Rd <- (cos(s/2))^2 + 0.1434143*steven
	Rd
}

Steven4 <- function(data)
{
	theta <- d2r(data$theta)
	z <- d2r(data$Z)
	s <- d2r(data$s)
	rb <- pmax(0, cos(theta)/cos(z))
	steven <- sin(s) - s*cos(s) -pi*(sin(s/2))^2
	#select a and b values
	ab.range <- d2r(c(35,45,55,65))
	tmp <- abs(z - matrix(rep(ab.range, length(z)), ncol=4, byrow = TRUE))
	choice <- apply(tmp, 1, which.min)
	a <- c(0.63, 0.6, 0.53, 0.46); b <- c(-1.04, -1, -0.9, -0.85);
	Rd <- ifelse(data$epsilon >= 1.1, a[choice]*rb + (1-a[choice])*((cos(s/2))^2 + 2*b[choice]/pi/(3 + 2*b[choice])*steven), (cos(s/2))^2 + 0.1434143*steven)
	Rd
}


#################################################################################
# Willmott 1982 (s, theta, z, I, Isc)
#################################################################################
Willmott <- function(data)
{
	AI <- data$I/1362
	theta <- d2r(data$theta)
	z <- d2r(data$Z)
	s <- d2r(data$s)
	rb <- pmax(0, cos(theta)/cos(z))
	cb <- pmin(pmax(1.0115-0.20293*s-0.080823*s^2, 0.5),1)
	Rd <- AI*rb + (1-AI)*cb
	Rd
}

#################################################################################
# Skartveit 1986 (s, theta, z, Ih, Ioh)
#################################################################################
Skartveit <- function(data)
{
	AI <- (data$Gh - data$Dh)/data$Ioh
	theta <- d2r(data$theta)
	z <- d2r(data$Z)
	s <- d2r(data$s)
	rb <- pmax(0, cos(theta)/cos(z))
	Zz <- pmax(0, 0.3-2*AI)
	Rd <- AI*rb + (1-AI-Zz)*(cos(s/2))^2 + Zz*cos(s)
	Rd
}

#################################################################################
# Koronakis 1986 (s)
#################################################################################
Koronakis <- function(data)
{
	s <- d2r(data$s)
	Rd <- (2+cos(s))/3
	Rd
}

#################################################################################
# Gueymard 1987 (s, z, theta, Gh, Dh)
#################################################################################
Gueymard <- function(data)
{
	theta <- d2r(data$theta)
	z <- d2r(data$Z)
	s <- d2r(data$s)
	h = 0.01*(90-data$Z)
	THETA <- data.matrix(data.frame(1, cos(theta), (cos(theta))^2, (cos(theta))^3))
	H <- t(data.matrix(data.frame(1, h, h^2, h^3, h^4)))
	X <- matrix(c(-0.897, -3.364, 3.960, -1.909, 0, 4.448, -12.962,  34.601, -48.784, 27.511, -2.770, 9.164, -18.876, 23.776, -13.014, 0.312, -0.217, -0.805, 0.318, 0), nrow = 4, byrow = TRUE)
	tmp <- rowSums(t(X%*%H)*THETA)
	Ff <- (1-0.2249*(sin(s))^2 + 0.1231*sin(2*s) - 0.0342*sin(4*s))/(1-0.2249)
	Gg <- 0.408-0.323*h+0.384*h^2 -0.17*h^3
	Rd0 <- exp(tmp) + Ff*Gg
	Y <- ifelse(data$Dh/data$Gh <= 0.227, 6.6667*data$Dh/data$Gh-1.4167, 1.2121*data$Dh/data$Gh-0.1758)
	Npt <- pmax(pmin(Y,1),0)
	b <- 0.5 +Npt
	Rd1 <- (cos(s/2))^2 + 2*b/(3+2*b)/pi*(sin(s) - s*cos(s) -pi*(sin(s/2))^2)
	Rd <- (1-Npt)*Rd0 + Npt*Rd1
	Rd
}

#################################################################################
# Muneer 1990 (s, theta, z, Gh, Dh, Ih, Ioh)
#################################################################################
Muneer1 <- function(data)
{
	AI <- (data$Gh - data$Dh)/data$Ioh
	theta <- d2r(data$theta)
	z <- d2r(data$Z)
	s <- d2r(data$s)
	rb <- pmax(0, cos(theta)/cos(z))
	steven <- sin(s) - s*cos(s) -pi*(sin(s/2))^2
	Rd <- ifelse(theta > 1.56207, (cos(s/2))^2+0.2522705*steven, ifelse(abs(data$Gh-data$Dh)<0.05*data$Gh, (cos(s/2))^2+ 0.1681637*steven, AI*rb + (1-AI)*((cos(s/2))^2-0.2242638*steven)))
	Rd
}

Muneer2 <- function(data)
{
	AI <- (data$Gh - data$Dh)/data$Ioh
	theta <- d2r(data$theta)
	z <- d2r(data$Z)
	s <- d2r(data$s)
	rb <- pmax(0, cos(theta)/cos(z))
	steven <- sin(s) - s*cos(s) -pi*(sin(s/2))^2
	Rd <- ifelse(theta > 1.56207, (cos(s/2))^2+0.2522705*steven, ifelse(abs(data$Gh-data$Dh)<0.05*data$Gh, (cos(s/2))^2+ 0.1681637*steven, AI*rb + (1-AI)*((cos(s/2))^2+(0.04-0.820*AI-2.0260*AI^2)*steven)))
	Rd
}

#################################################################################
# Reindl 1990 (s, z, theta, Ih, Ioh, Gh)
#################################################################################
Reindl <- function(data)
{
	AI <- (data$Gh - data$Dh)/data$Ioh
	theta <- d2r(data$theta)
	z <- d2r(data$Z)
	s <- d2r(data$s)
	rb <- pmax(0, cos(theta)/cos(z))
	ff <- sqrt((data$Gh - data$Dh)/data$Gh)
	Rd <- (1-AI)*(cos(s/2))^2*(1+ff*(sin(s/2))^3) + AI*rb
	Rd
}

#################################################################################
# Perez 1990 (s, z, theta, Ih, Ioh, Gh)
#################################################################################
Perez1 <- function(data) #Eq.(8) Perez 1987
{
  Gh = data$Gh
  Dh = data$Dh
  Ih = Gh-Dh
  theta <- d2r(data$theta)
  z <- d2r(data$Z)
  s <- d2r(data$s)
  Delta = data$Delta
  epsilon <- pmin((Ih/cos(z)+Dh)/Dh, 999) #note the epsilon here is different from Perez1990
  alpha <- d2r(25)
  xi <- d2r(6.5)

  #a, b, c, d
  psi.h <- ifelse(z > pi/2 - alpha, (pi/2-z + alpha)/2/alpha, 1)
  psi.c <- (pi/2-theta+alpha)/2/alpha
  chi.c <- ifelse(theta < pi/2-alpha, psi.h*cos(theta), ifelse(theta >pi/2-alpha & theta < pi/2, psi.h*psi.c*sin(psi.c*alpha), 0))
  a <- 2*(1-cos(alpha))*chi.c

  chi.h <- ifelse(z < pi/2 - alpha, cos(z), psi.h*sin(psi.h*alpha))
  c <- 2*(1-cos(alpha))*chi.h

  #assign points to the correct epsilon bin
  e.bin <- load.Perez.par(model = "Perez87a")$epsilon.mat #Perez3 shares the same epsilon range as Perez2
  tmp <- lapply(1:8, function(i) which(epsilon >= e.bin[i,1] & epsilon < e.bin[i,2]))
  choice <- array(NA, length(z))
  for(i in 1:8){choice[tmp[[i]]] = i}
  Par <- load.Perez.par(model = "Perez87a")$F.mat[choice,]

  F1 <- pmax(Par[,1]+Delta*Par[,2]+z*Par[,3], 0)
  F2 <- (Par[,4]+Delta*Par[,5]+z*Par[,6])
  Rd = (1-F1)*0.5*(1+cos(s)) + F1*(a/c) + F2*sin(s)
  Rd
}

Perez2 <- function(data)
{
  Gh = data$Gh
  Dh = data$Dh
  Ih = Gh-Dh
	theta <- d2r(data$theta)
	z <- d2r(data$Z)
	s <- d2r(data$s)
	Delta = data$Delta
	rb <- pmax(0, cos(theta)/cos(z))
	epsilon <- pmin((Ih/cos(z)+Dh)/Dh, 999) #note the epsilon here is different from Perez1990

	#assign points to the correct epsilon bin
	e.bin <- load.Perez.par(model = "Perez87b")$epsilon.mat
	tmp <- lapply(1:8, function(i) which(epsilon >= e.bin[i,1] & epsilon < e.bin[i,2]))
	choice <- array(NA, length(z))
	for(i in 1:8){choice[tmp[[i]]] = i}
	Par <- load.Perez.par(model = "Perez87b")$F.mat[choice,]

	F1 <- pmax(Par[,1]+Delta*Par[,2]+z*Par[,3], 0)
	F2 <- (Par[,4]+Delta*Par[,5]+z*Par[,6])
	Rd = (1-F1)*0.5*(1+cos(s)) + F1*rb + F2*sin(s)
	Rd
}

Perez3 <- function(data)
{
  theta <- d2r(data$theta)
  z <- d2r(data$Z)
  s <- d2r(data$s)
  Delta = data$Delta
  rb <- pmax(0, cos(theta)/cos(z))
  epsilon <- pmin(data$epsilon, 11.999)

  #assign points to the correct epsilon bin
  e.bin <- load.Perez.par(model = "Perez90")$epsilon.mat
  tmp <- lapply(1:8, function(i) which(epsilon >= e.bin[i,1] & epsilon < e.bin[i,2]))
  choice <- array(NA, length(z))
  for(i in 1:8){choice[tmp[[i]]] = i}
  Par <- load.Perez.par(model = "Perez90")$F.mat[choice,]

  F1 <- pmax(Par[,1]+Delta*Par[,2]+z*Par[,3], 0)
  F2 <- (Par[,4]+Delta*Par[,5]+z*Par[,6])
  Rd = (1-F1)*0.5*(1+cos(s)) + F1*rb + F2*sin(s)
  Rd
}

Perez4 <- function(data) #Eq.(8) Perez 1987
{
	theta <- d2r(data$theta)
	z <- d2r(data$Z)
	s <- d2r(data$s)
	Delta = data$Delta
	epsilon <- pmin(data$epsilon, 11.999)
	alpha <- d2r(25)
	xi <- d2r(6.5)

	#a, b, c, d
	psi.h <- ifelse(z > pi/2 - alpha, (pi/2-z + alpha)/2/alpha, 1)
	psi.c <- (pi/2-theta+alpha)/2/alpha
	chi.c <- ifelse(theta < pi/2-alpha, psi.h*cos(theta), ifelse(theta >pi/2-alpha & theta < pi/2, psi.h*psi.c*sin(psi.c*alpha), 0))
	a <- 2*(1-cos(alpha))*chi.c

	chi.h <- ifelse(z < pi/2 - alpha, cos(z), psi.h*sin(psi.h*alpha))
	c <- 2*(1-cos(alpha))*chi.h

	#assign points to the correct epsilon bin
	e.bin <- load.Perez.par(model = "Perez90")$epsilon.mat #Perez3 shares the same epsilon range as Perez2
	tmp <- lapply(1:8, function(i) which(epsilon >= e.bin[i,1] & epsilon < e.bin[i,2]))
	choice <- array(NA, length(z))
	for(i in 1:8){choice[tmp[[i]]] = i}
	Par <- load.Perez.par(model = "Perez90")$F.mat[choice,]

	F1 <- pmax(Par[,1]+Delta*Par[,2]+z*Par[,3], 0)
	F2 <- (Par[,4]+Delta*Par[,5]+z*Par[,6])
	Rd = (1-F1)*0.5*(1+cos(s)) + F1*(a/c) + F2*sin(s)
	Rd
}

#################################################################################
# Olmo et al. 1999 (s, z, theta, Ih, Ioh, Gh)
#################################################################################
Olmo1 <- function(data)
{
  Gh = data$Gh
  Dh = data$Dh
  Ih = Gh-Dh
  theta <- d2r(data$theta)
  z <- d2r(data$Z)
  s <- d2r(data$s)
  rho = data$rho
  rb <- pmax(0, cos(theta)/cos(z))
  #Ic <- cos(d2r(data$theta))/cos(d2r(data$Z))*(data$Gh-data$Dh)
  #Dg <- data$Gh*0.5*data$rho*(1-cos(d2r(data$s)))

  fc <- 1+rho*(sin(theta/2))^2
  kt <- Gh/data$Ioh
  Rd = Gh*(fc*exp(-kt*(theta^2-z^2)) - 0.5*data$rho*(1-cos(s)))/Dh - Ih*rb/Dh
  #Rd = Gh*(fc*exp(-kt*(theta^2-z^2)))/Dh - Dg/Dh - Ic/Dh
  Rd
}

Olmo2 <- function(data)
{
  Gh = data$Gh
  theta <- d2r(data$theta)
  z <- d2r(data$Z)
  rho = data$rho
  kt <- Gh/data$Ioh

  fc <- ifelse(kt > 0.65, 1, ifelse(kt < 0.35, 1-rho*(cos(theta/2))^3, 1-rho*sin(theta/2)))
  Rd = fc*exp(-kt*(theta^2-z^2))
  Rd
}

#################################################################################
# Tian 2001 (s)
#################################################################################
Tian <- function(data)
{
	Rd <- 1-data$s/180
	Rd
}

#################################################################################
# Badescu 2002 (s)
#################################################################################
Badescu <- function(data)
{
	s <- d2r(data$s)
	Rd <- (3+cos(2*s))/4
	Rd
}

#################################################################################
# Error calculation functions
#################################################################################
MBE <- function(meas,pred)
{
	mean(pred-meas)/mean(meas)*100
}
RMSE <- function(meas,pred)
{
	sqrt(mean((pred-meas)^2))/mean(meas)*100
}

