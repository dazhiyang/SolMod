# Transposition modeling of solar irradiance

Transposition models convert solar irradiance received by a horizontal surface to that received by an arbitrary tilted surface. The inverse transposition models are used to convert irradiance from tilt to horizontal. This package gives several popular models to perform such bidirectional conversion. 

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

### Prerequisites

This is an R package, so you need to install [R](https://www.r-project.org/) on your computer first. In addition, [RStudio](https://www.rstudio.com/) is an integrated development environment (IDE) for R; it is highly recommended.

### Installing

Once R and RStudio are installed. Open R or RStudio and install the [devtools](https://cran.r-project.org/web/packages/devtools/index.html) package, which allows you to install R package from GitHub

```
install.packages("devtools")
```

Load the package that you just installed

```
library("devtools")
```

Now, you can install the SolMod package, using

```
install_github("dazhiyang/SolMod")
```

## Running the tests

This code segment gives an example on how to run transposition modeling (horizontal to tilt) using a variety of models. 

```
data("NREL") #load data
slope = c(40, 90, 90, 90, 90) # initialize tilt angles
azimuth.PV = c(180, 0,  90, 180, 270) # initialize azimuth angles

Gc <- list()
for(i in 1:length(slope)) #loop for the 5 tilts
{
	NREL$theta <- NREL[,which(substr(names(NREL), 1, 2)=="th")[i]]
	NREL$s <- slope[i]

  	Rd <- data.frame(tm = NREL$Tm) #initialize a data.frame

	#there are a total of 26 models, use 3 for demonstration purposes
  	Rd$Liu <- Liu(NREL) #Liu
      	Rd$Perez4 <- Perez4(NREL) #Perez4
  	Rd$Gueymard <- Gueymard(NREL) #Gueymard

	#compute tilted irradiance
  	Rd <- Rd[,-1]
  	Ic <- cos(d2r(NREL$theta))/cos(d2r(NREL$Z))*(NREL$Gh-NREL$Dh)
  	Dg <- NREL$Gh*0.5*NREL$rho*(1-cos(d2r(NREL$s)))
  	Dc <- Rd*NREL$Dh
  	Gc[[i]] <- Ic + Dg + Dc
}

#save predicted and measured Gc values
pred <- Gc
meas <- NREL[,which(substr(names(NREL), 1, 2)=="Gc")]

#error calculation
rmse <- matrix(unlist(lapply(1:length(slope), function(y) apply(Gc[[y]], 2, function(x) round(RMSE(x, meas[,y]),1)))), ncol = length(slope))
rownames(rmse) <- names(Gc[[1]])
colnames(rmse) <- paste(" (", slope, ", ", azimuth.PV, ")", sep = "")
rmse
```

## For more information

* The original paper, a review published in Solar Energy, describing these transposition models can be found [here](https://doi.org/10.1016/j.solener.2016.06.062).
* The original document describing the best performing model, namely, the Perez model, can be found [here](https://doi.org/10.1016/0038-092X(90)90055-H).


## License

This package is free and open source software, licensed under GPL-3.
