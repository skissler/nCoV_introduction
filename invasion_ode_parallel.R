# =============================================================================
# Import
# =============================================================================

library(tidyverse)
library(deSolve)
library(lubridate)
library(parallel)
library(lhs)

nrevssCDC_ILI <- read.csv('nrevssCDC_ILI.csv') %>% as_tibble()
nrevssCDC_ILI$WEEKEND <- as_date(ymd(nrevssCDC_ILI$WEEKEND))

ssecalc <- function(scalingfactor, odefit){
	odefit %>% 
		filter(time>=extractrange[1] & time<=extractrange[2]) %>%
		mutate(Strain1=I1S2 + I1E2 + I1I2 + I1R2) %>%
		mutate(Strain2=S1I2 + E1I2 + I1I2 + R1I2) %>%
		mutate(WEEKINDEX=1:nrow(.)) %>%
		left_join(nrevssCDC_ILI, by="WEEKINDEX") %>%
		mutate(SSE=(scalingfactor*Strain1-OC43_ILI)^2 + (scalingfactor*Strain2-HKU1_ILI)^2) %>%
		summarise(SSE=sum(SSE)) %>%
		as.numeric()
}

# =============================================================================
# Define and simulate model
# =============================================================================

source('invasion_rates.R')

model <- function(t, y, parms) {
		dS1S2 <- -S1S2c1(t,y,parms) - S1S2c2(t,y,parms) + 
				  R1S2c1(t,y,parms) + S1R2c2(t,y,parms)
		dE1S2 <- -E1S2c1(t,y,parms) - E1S2c2(t,y,parms) + 
				  S1S2c1(t,y,parms) + E1R2c2(t,y,parms)
		dS1E2 <- -S1E2c1(t,y,parms) - S1E2c2(t,y,parms) + 
				  R1E2c1(t,y,parms) + S1S2c2(t,y,parms)
		dE1E2 <- -E1E2c1(t,y,parms) - E1E2c2(t,y,parms) + 
				  S1E2c1(t,y,parms) + E1S2c2(t,y,parms)
		dI1S2 <- -I1S2c1(t,y,parms) - I1S2c2(t,y,parms) + 
				  E1S2c1(t,y,parms) + I1R2c2(t,y,parms)
		dS1I2 <- -S1I2c1(t,y,parms) - S1I2c2(t,y,parms) + 
				  R1I2c1(t,y,parms) + S1E2c2(t,y,parms)
		dR1S2 <- -R1S2c1(t,y,parms) - R1S2c2(t,y,parms) + 
				  I1S2c1(t,y,parms) + R1R2c2(t,y,parms)
		dI1E2 <- -I1E2c1(t,y,parms) - I1E2c2(t,y,parms) + 
				  E1E2c1(t,y,parms) + I1S2c2(t,y,parms)
		dE1I2 <- -E1I2c1(t,y,parms) - E1I2c2(t,y,parms) + 
				  S1I2c1(t,y,parms) + E1E2c2(t,y,parms)
		dS1R2 <- -S1R2c1(t,y,parms) - S1R2c2(t,y,parms) + 
				  R1R2c1(t,y,parms) + S1I2c2(t,y,parms)
		dR1E2 <- -R1E2c1(t,y,parms) - R1E2c2(t,y,parms) + 
				  I1E2c1(t,y,parms) + R1S2c2(t,y,parms)
		dI1I2 <- -I1I2c1(t,y,parms) - I1I2c2(t,y,parms) + 
				  E1I2c1(t,y,parms) + I1E2c2(t,y,parms)
		dE1R2 <- -E1R2c1(t,y,parms) - E1R2c2(t,y,parms) + 
				  S1R2c1(t,y,parms) + E1I2c2(t,y,parms)
		dR1I2 <- -R1I2c1(t,y,parms) - R1I2c2(t,y,parms) + 
				  I1I2c1(t,y,parms) + R1E2c2(t,y,parms)
		dI1R2 <- -I1R2c1(t,y,parms) - I1R2c2(t,y,parms) + 
				  E1R2c1(t,y,parms) + I1I2c2(t,y,parms)
		dR1R2 <- -R1R2c1(t,y,parms) - R1R2c2(t,y,parms) + 
				  I1R2c1(t,y,parms) + R1I2c2(t,y,parms)
		return(list(c(dS1S2,dE1S2,dS1E2,dE1E2,dI1S2,dS1I2,dR1S2,dI1E2,dE1I2,dS1R2,dR1E2,dI1I2,dE1R2,dR1I2,dI1R2,dR1R2)))
}

y <- c(S1S2 = 1,
	   E1S2 = 0,
	   S1E2 = 0,
	   E1E2 = 0,
	   I1S2 = 0,
	   S1I2 = 0,
	   R1S2 = 0,
	   I1E2 = 0,
	   E1I2 = 0,
	   S1R2 = 0,
	   R1E2 = 0,
	   I1I2 = 0,
	   E1R2 = 0,
	   R1I2 = 0,
	   I1R2 = 0,
	   R1R2 = 0) # Initial conditions

times <- seq(0,52*30,1) # Time range for simulation, weeks

extractrange <- c(52*24.5, 52*29.5) # Which weeks to compare with data? 

# Set test values for the parameters
variable.parmranges <- list("sigma1.val"=c(1/100, 1/25), 
						 "sigma2.val"=c(1/100, 1/25), 
						 "nu.val"=c(0.5, 2),
						 "gamma.val"=c(0.5, 2),
						 "chi12.val"=c(0, 1), 
						 "chi21.val"=c(0, 1), 
						 "amplitude"=c(0, 1), 
						 "baseline"=c(1, 2), 
						 "phi.val"=c(-8, 8))

constant.parmvals <- list(
						  "kappa.val"=0.01,
						  "importtime1"=52,
						  "importtime2"=0,
						  "importlength"=0.5)

ndraws <- 100000

parmsdf.var <- as.data.frame(randomLHS(ndraws, length(variable.parmranges))) 
names(parmsdf.var) <- names(variable.parmranges)
for(parm in names(variable.parmranges)){
	parmsdf.var[,parm] <- parmsdf.var[,parm]*(variable.parmranges[[parm]][2]-variable.parmranges[[parm]][1]) + variable.parmranges[[parm]][1]
}

parmsdf.const <- data.frame(init=1:ndraws)
for(parm in names(constant.parmvals)){
	parmsdf.const[,parm] <- rep(constant.parmvals[[parm]], ndraws)
}
parmsdf.const$init <- NULL

parmsdf <- cbind(parmsdf.var, parmsdf.const)
rm(parmsdf.var)
rm(parmsdf.const)

# Update the import times (comment if you want to keep them constant): 
parmsdf$importtime1 <- rep(c(52,52,0), ndraws)[1:ndraws]
parmsdf$importtime2 <- rep(c(52,0,52), ndraws)[1:ndraws]

parmslist <- lapply(split(parmsdf, seq(nrow(parmsdf))),unlist)

# -----------------------------------------------------------------------------

scalingfactors <- seq(0.001, 0.1, 0.001) 

getsse <- function(parms, y, times, model){

	odefit <- tryCatch({as_tibble(as.data.frame(lsoda(y,times,model,parms)))}, 
	warning=function(w){print("A warning occurred"); data.frame(time=NA, output=NA)},
	error=function(e){print("An error occurred"); data.frame(time=NA, output=NA)})

	if(nrow(odefit)<length(times)){
		SSEvals <- data.frame(SF=NA, SSE=NA)
		return(SSEvals)
		}else{
		SSEvals <- data.frame(SF=scalingfactors, SSE=unlist(lapply(scalingfactors, ssecalc, odefit))) %>% 
			filter(SSE==min(.$SSE))
		return(SSEvals[1,])
	}
}


# Run the simulations in parallel: 
start_time = Sys.time()
SSEdf <- bind_rows(mclapply(parmslist, getsse, y, times, model, mc.cores=20))
end_time = Sys.time()

# Add a column with the SSEs: 
parmsdf <- cbind(parmsdf, SSEdf)

# Save the output: 
write.csv(parmsdf, file="parmsdf_LHS_100000.csv", row.names=FALSE)

