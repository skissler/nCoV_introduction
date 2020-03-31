# =============================================================================
# Import
# =============================================================================

library(tidyverse) 
library(deSolve)
library(lubridate)

# nrevss_beta <- read.csv('data/nrevss_beta.csv') %>% as_tibble()
# nrevss_beta$WEEKEND <- as_date(ymd(nrevss_beta$WEEKEND))

nrevssCDC_ILI <- read.csv('data/nrevssCDC_ILI.csv') %>% as_tibble()
nrevssCDC_ILI$WEEKEND <- as_date(ymd(nrevssCDC_ILI$WEEKEND))

# =============================================================================
# Define and simulate model
# =============================================================================

source('code/invasion_rates.R')

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

parms=c(
	sigma1.val = 1/40, # Waning immunity rate, strain 1, weeks (def. 1/40)
	sigma2.val = 1/40, # Waning immunity rate, strain 2, weeks (def. 1/40)
	nu.val = 1/(5/7),      # Rate of progression to infection, weeks (def. 1/1)
	gamma.val = 1/(5/7),   # Rate of recovery, weeks (def. 1/1)
	chi12.val = 0.7,   # Cross immunity, strain 1 against strain 2 (def. 0.7)
	chi21.val = 0.5,   # Cross immunity, strain 2 against strain 1 (def. 0.3)
	amplitude = 0.6,     # R0 seasonality amplitude (dev. 1)
	baseline = 1.4,    # R0 seasonality baseline (def. 1.5)
	phi.val =  -4,      # R0 seasonality phase shift, weeks (def. -4)
	kappa.val = 0.01,  # Force of infection pulse size for introductions  
	importtime1 = 0,  # Time of strain 1 importation, weeks
	importtime2 = 52,   # Time of strain 2 importation, weeks
	importlength = 0.5 # Duration of importation pulse, weeks
	)

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
	   R1R2 = 0)

times <- seq(0,52*30,1)

odefit <- as.data.frame(lsoda(y,times,model,parms)) %>% as_tibble()

# =============================================================================
# Calculate SSE
# =============================================================================

extractrange <- c(52*24.5, 52*29.5)
scalingfactors <- seq(0.001, 0.1, 0.001)

getSSE <- function(scalingfactor, odefit){
	odefit %>% 
		filter(time>=extractrange[1] & time<=extractrange[2]) %>%
		mutate(Strain1=I1S2 + I1E2 + I1I2 + I1R2) %>%
		mutate(Strain2=S1I2 + E1I2 + I1I2 + R1I2) %>%
		mutate(WEEKINDEX=1:n()) %>%
		left_join(nrevssCDC_ILI, by="WEEKINDEX") %>%
		mutate(SSE=(scalingfactor*Strain1-OC43_ILI)^2 + (scalingfactor*Strain2-HKU1_ILI)^2) %>%
		summarise(SSE=sum(SSE)) %>%
		as.numeric()
}

SSEvals <- data.frame(SF=scalingfactors, SSE=unlist(lapply(scalingfactors, getSSE, odefit)))

SSE <- SSEvals %>% 
	filter(SSE==min(.$SSE))

# SSE <- odefit %>% 
# 	filter(time>=extractrange[1] & time<=extractrange[2]) %>%
# 	mutate(Strain1=I1S2 + I1E2 + I1I2 + I1R2) %>%
# 	mutate(Strain2=S1I2 + E1I2 + I1I2 + R1I2) %>%
# 	mutate(WEEKINDEX=1:nrow(.)) %>%
# 	left_join(nrevssCDC_ILI, by="WEEKINDEX") %>%
# 	mutate(SSE=(scalingfactor*Strain1-OC43_ILI)^2 + (scalingfactor*Strain2-HKU1_ILI)^2) %>%
# 	summarise(SSE=sum(SSE)) %>%
# 	as.numeric()

# =============================================================================
# Plot the simulation in a custom window
# =============================================================================

# plotrange <- c(0, 52*30)
plotrange <- c(52*25, 52*30)

figodefit <- odefit %>% 
	filter(time>=plotrange[1] & time<=plotrange[2]) %>%
	mutate(Strain1=I1S2 + I1E2 + I1I2 + I1R2) %>%
	mutate(Strain2=S1I2 + E1I2 + I1I2 + R1I2) %>%
	mutate(time_yr=time/52) %>%
	select(time_yr, Strain1, Strain2) %>%
	pivot_longer(c("Strain1","Strain2"), names_to="Virus", values_to="ProportionInfected") %>%
	ggplot(aes(x=time_yr, y=ProportionInfected, col=Virus)) + 
		geom_line(size=0.8, alpha=0.6) + 
		scale_x_continuous(breaks=0:100) + 
		scale_color_manual(values=c("blue","red")) + 
		labs(x="Year", y="Proportion Infected") + 
		theme_minimal() + 
		theme(text=element_text(size=14))
# ggsave(figodefit, file="/Users/sk792/DropboxHarvard/Projects/Wuhan/figures/odefit.pdf", width=8, height=5)


# =============================================================================
# Plot susceptibles in a custom window
# =============================================================================

plotrange <- c(52*24.5, 52*29.5)

figodefit_susc <- odefit %>% 
	filter(time>=plotrange[1] & time<=plotrange[2]) %>%
	mutate(Strain1=(S1S2 + (1-parms[["chi21.val"]])*(S1E2 + S1I2 + S1R2))) %>%
	mutate(Strain2=(S1S2 + (1-parms[["chi12.val"]])*(E1S2 + I1S2 + R1S2))) %>%
	mutate(time_yr=time/52) %>%
	select(time_yr, Strain1, Strain2) %>%
	pivot_longer(c("Strain1","Strain2"), names_to="Virus", values_to="ProportionSusceptible") %>%
	ggplot(aes(x=time_yr, y=ProportionSusceptible, col=Virus)) + 
		geom_line(size=0.8, alpha=0.6) + 
		scale_x_continuous(breaks=0:100) + 
		scale_color_manual(values=c("blue","red")) + 
		labs(x="Year", y="Proportion Susceptible") + 
		theme_minimal() + 
		expand_limits(y=c(0,1)) + 
		theme(text=element_text(size=14))

# =============================================================================
# Plot simulation against NREVSS percent positive data
# =============================================================================

scalingfactor <- 1

figodefit_nrevss_pp <- odefit %>% 
	filter(time>=extractrange[1] & time<=extractrange[2]) %>%
	mutate(Strain1=I1S2 + I1E2 + I1I2 + I1R2) %>%
	mutate(Strain2=S1I2 + E1I2 + I1I2 + R1I2) %>%
	mutate(WEEKINDEX=1:n()) %>%
	left_join(nrevssCDC_ILI, by="WEEKINDEX") %>%
	select(WEEKEND, Strain1, Strain2, OC43, HKU1) %>%
	mutate(Strain1=Strain1*scalingfactor) %>% 
	mutate(Strain2=Strain2*scalingfactor) %>% 
	pivot_longer(c("Strain1","Strain2","OC43","HKU1"), names_to="Virus", values_to="PercentPositive") %>%
	mutate(DataType=case_when(Virus%in%c("OC43","HKU1")~"Actual",TRUE~"Simulated")) %>%
	ggplot(aes(x=WEEKEND, y=PercentPositive, col=Virus, linetype=DataType)) + 
		geom_line(size=0.8, alpha=0.6) + 
		scale_color_manual(values=c("red","blue","blue","red")) + 
		labs(x="Year", y="Proportion Infected") + 
		theme_minimal() + 
		theme(text=element_text(size=14))
# ggsave(figodefit_nrevss_pp, file="/Users/sk792/DropboxHarvard/Projects/Wuhan/figures/odefit_nrevss_pp.pdf", width=8, height=5)

# =============================================================================
# Plot simulation against NREVSS PP*ILI data
# =============================================================================

scalingfactor <- 0.075

figodefit_nrevss <- odefit %>% 
	filter(time>=extractrange[1] & time<=extractrange[2]) %>%
	mutate(Strain1=I1S2 + I1E2 + I1I2 + I1R2) %>%
	mutate(Strain2=S1I2 + E1I2 + I1I2 + R1I2) %>%
	mutate(WEEKINDEX=1:n()) %>%
	left_join(nrevssCDC_ILI, by="WEEKINDEX") %>%
	select(WEEKEND, Strain1, Strain2, OC43_ILI, HKU1_ILI) %>%
	mutate(OC43_ILI=OC43_ILI) %>%
	mutate(HKU1_ILI=HKU1_ILI) %>% 
	mutate(Strain1=Strain1*scalingfactor) %>% 
	mutate(Strain2=Strain2*scalingfactor) %>% 
	pivot_longer(c("Strain1","Strain2","OC43_ILI","HKU1_ILI"), names_to="Virus", values_to="ProportionInfected") %>%
	mutate(DataType=case_when(Virus%in%c("OC43_ILI","HKU1_ILI")~"Actual",TRUE~"Simulated")) %>%
	ggplot(aes(x=WEEKEND, y=ProportionInfected, col=Virus, linetype=DataType)) + 
		geom_line(size=0.8, alpha=0.6) + 
		scale_color_manual(values=c("red","blue","blue","red")) + 
		labs(x="Year", y="Proportion Infected") + 
		theme_minimal() + 
		theme(text=element_text(size=14))
# ggsave(figodefit_nrevss, file="/Users/sk792/DropboxHarvard/Projects/Wuhan/figures/odefit_nrevss_div20.pdf", width=8, height=5)

# =============================================================================
# Plot effective reproduction number
# =============================================================================

plotrange <- c(52*24.5, 52*29.5)

figRe <- odefit %>% 
	filter(time>=plotrange[1] & time<=plotrange[2]) %>%
	mutate(Strain1=(S1S2 + (1-parms[["chi21.val"]])*(S1E2 + S1I2 + S1R2))*beta.val(time, parms[["amplitude"]], parms[["baseline"]], parms[["phi.val"]])) %>%
	mutate(Strain2=(S1S2 + (1-parms[["chi12.val"]])*(E1S2 + I1S2 + R1S2))*beta.val(time, parms[["amplitude"]], parms[["baseline"]], parms[["phi.val"]])) %>%
	mutate(time_yr=time/52+2015-25) %>%
	select(time_yr, Strain1, Strain2) %>%
	pivot_longer(c("Strain1","Strain2"), names_to="Virus", values_to="Re") %>%
	ggplot(aes(x=time_yr, y=Re, col=Virus)) + 
		geom_line(size=0.8, alpha=0.6) + 
		scale_x_continuous(breaks=0:10000) + 
		scale_color_manual(values=c("blue","red")) + 
		labs(x="Year", y="Re") + 
		theme_minimal() + 
		theme(text=element_text(size=14))
# ggsave(figRe, file="/Users/sk792/DropboxHarvard/Projects/Wuhan/figures/Re.pdf", width=8, height=5)

# =============================================================================
# Plot effective reproduction number against Christine's data
# =============================================================================

wtvals <- read.csv("/Users/sk792/DropboxHarvard/Projects/Wuhan/data/2020.02.11R_ili_x_pos_pct_SARS.csv") %>% as_tibble()
wtvals$Week_start <- ymd(wtvals$Week_start)
wtvals$WEEKEND = wtvals$Week_start + days(6)
wtvals <- wtvals %>% 
	select(WEEKEND, CoVHKU1_R_ili_x_pos_pct_SARS, CoVOC43_R_ili_x_pos_pct_SARS) %>%
	rename(Re_OC43 = CoVOC43_R_ili_x_pos_pct_SARS) %>%
	rename(Re_HKU1 = CoVHKU1_R_ili_x_pos_pct_SARS)

wtvals %>% 
	pivot_longer(-WEEKEND) %>%
	ggplot(aes(x=WEEKEND, y=value, col=name)) + 
		geom_point(alpha=0.3) + 
		geom_line(stat="smooth", method="loess", span=0.2, size=1.2) + 
		theme_minimal()

redf_sim <- odefit %>% 
	filter(time>=extractrange[1] & time<=extractrange[2]) %>%
	mutate(Strain1=(S1S2 + (1-parms[["chi21.val"]])*(S1E2 + S1I2 + S1R2))*beta.val(time, parms[["amplitude"]], parms[["baseline"]], parms[["phi.val"]])) %>%
	mutate(Strain2=(S1S2 + (1-parms[["chi12.val"]])*(E1S2 + I1S2 + R1S2))*beta.val(time, parms[["amplitude"]], parms[["baseline"]], parms[["phi.val"]])) %>%
	mutate(WEEKINDEX=1:nrow(.)) %>%
	left_join(nrevssCDC_ILI, by="WEEKINDEX") %>%
	select(WEEKEND, Strain1, Strain2) %>%
	rename(OC43=Strain1, HKU1=Strain2) %>%
	pivot_longer(-WEEKEND, names_to="Virus", values_to="Re") 

redf_wt <- pivot_longer(wtvals,-WEEKEND, names_to="Virus", values_to="Re") %>%
	mutate(Virus=case_when(Virus=="Re_HKU1"~"HKU1", Virus=="Re_OC43"~"OC43"))

figwtre <- ggplot() + 
	geom_line(data=redf_sim, aes(x=WEEKEND, y=Re, col=Virus), size=1.2, alpha=0.6) + 
	geom_point(data=redf_wt, aes(x=WEEKEND, y=Re, col=Virus), alpha=0.4) +
	# geom_line(data=redf_wt, aes(x=WEEKEND, y=Re, col=Virus), stat="smooth", method="loess", span=0.2, alpha=0.7, size=1) +
	scale_color_manual(limits=c("OC43","HKU1"), values=c("OC43"="blue", "HKU1"="orange")) + 
	labs(x="Week", y="Re") + 
	theme_minimal() + 
	theme(legend.position="none", text=element_text(size=20)) + 
	facet_wrap(~Virus, nrow=1)
# ggsave(figwtre, file="/Users/sk792/DropboxHarvard/Projects/Wuhan/writeup/seasonality/figures/wtre.pdf", width=12, height=5)
# ggsave(figwtre, file="/Users/sk792/DropboxHarvard/Projects/Wuhan/writeup/seasonality/figures/wtre.png", width=12, height=5)

