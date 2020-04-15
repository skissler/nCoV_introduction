# =============================================================================
# Copyright
# =============================================================================
# COPYRIGHT 2020 STEPHEN M. KISSLER
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.

#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.

#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.

# =============================================================================
# Import
# =============================================================================

library(tidyverse) 
library(deSolve)
library(lubridate) 

nrevssCDC_ILI <- read.csv('data/nrevssCDC_ILI.csv') %>% as_tibble()
nrevssCDC_ILI$WEEKEND <- as_date(ymd(nrevssCDC_ILI$WEEKEND))

# =============================================================================
# Plot the four strains at once: 
# =============================================================================

figfourstrain <- nrevssCDC_ILI %>% 
	select(WEEKEND, OC43_ILI, HKU1_ILI, NL63_ILI, E229_ILI) %>% 
	rename("OC43"=OC43_ILI, "HKU1"=HKU1_ILI, "NL63"=NL63_ILI, "229E"=E229_ILI) %>%
	pivot_longer(-WEEKEND, names_to="Virus", values_to="PPxPILI") %>% 
	mutate(Genus=case_when(Virus %in% c("NL63","229E")~"Alpha", TRUE~"Beta")) %>% 
	ggplot(aes(x=WEEKEND, y=PPxPILI, col=Virus)) + 
		geom_line(size=1.2, alpha=0.6) +
		scale_color_manual(limits=c("NL63","229E","OC43","HKU1"), values=c("OC43"="blue", "HKU1"="orange", "NL63"="red", "229E"="black")) + 
		labs(x="Week", y="% positive x % ILI") + 
		theme_minimal() + 
		theme(text=element_text(size=24)) + 
		facet_wrap(~Genus, nrow=1)
# ggsave(figfourstrain, file="/Users/sk792/DropboxHarvard/Projects/Wuhan/writeup/seasonality/figures/fourstrain.pdf", width=12, height=5)
# ggsave(figfourstrain, file="/Users/sk792/DropboxHarvard/Projects/Wuhan/writeup/seasonality/figures/fourstrain.png", width=12, height=5)


figalphabeta <- nrevssCDC_ILI %>% 
	select(WEEKEND, OC43_ILI, HKU1_ILI, NL63_ILI, E229_ILI) %>% 
	rename("OC43"=OC43_ILI, "HKU1"=HKU1_ILI, "NL63"=NL63_ILI, "229E"=E229_ILI) %>%
	mutate(Alpha = OC43 + HKU1, Beta = NL63 + `229E`) %>%
	select(WEEKEND, Alpha, Beta) %>%
	pivot_longer(-WEEKEND, names_to="Genus", values_to="PPxPILI") %>%
	ggplot(aes(x=WEEKEND, y=PPxPILI, col=Genus)) + 
		geom_line(size=1.2, alpha=0.6) +
		scale_color_manual(values=c("Alpha"="black", "Beta"="blue")) + 
		labs(x="Week", y="% positive x % ILI") + 
		theme_minimal() + 
		theme(text=element_text(size=24)) 
# ggsave(figalphabeta, file="/Users/sk792/DropboxHarvard/Projects/Wuhan/writeup/seasonality/figures/alphabeta.pdf", width=8, height=5)
# ggsave(figalphabeta, file="/Users/sk792/DropboxHarvard/Projects/Wuhan/writeup/seasonality/figures/alphabeta.png", width=8, height=5)

# =============================================================================
# Plot a good fit to the data and the associated Re:
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

# parms=c(
# 	sigma1.val = 1/40, # Waning immunity rate, strain 1, weeks (def. 1/40)
# 	sigma2.val = 1/40, # Waning immunity rate, strain 2, weeks (def. 1/40)
# 	nu.val = 1/(5/7),      # Rate of progression to infection, weeks (def. 1/1)
# 	gamma.val = 1/(5/7),   # Rate of recovery, weeks (def. 1/1)
# 	chi12.val = 0.7,   # Cross immunity, strain 1 against strain 2 (def. 0.7)
# 	chi21.val = 0.5,   # Cross immunity, strain 2 against strain 1 (def. 0.3)
# 	amplitude = 0.6,     # R0 seasonality amplitude (dev. 1)
# 	baseline = 1.4,    # R0 seasonality baseline (def. 1.5)
# 	phi.val = -4,      # R0 seasonality phase shift, weeks (def. -4)
# 	kappa.val = 0.01,  # Force of infection pulse size for introductions  
# 	importtime1 = 0,  # Time of strain 1 importation, weeks
# 	importtime2 = 52,   # Time of strain 2 importation, weeks
# 	importlength = 0.5 # Duration of importation pulse, weeks
# 	)

# These are the ones from the original paper: 
# parms=c(
# 	sigma1.val = 1/40, # Waning immunity rate, strain 1, weeks (def. 1/40)
# 	sigma2.val = 1/38, # Waning immunity rate, strain 2, weeks (def. 1/40)
# 	nu.val = 1/(5/7),      # Rate of progression to infection, weeks (def. 1/1)
# 	gamma.val = 1/(4.9/7),   # Rate of recovery, weeks (def. 1/1)
# 	chi12.val = 0.74,   # Cross immunity, strain 1 against strain 2 (def. 0.7)
# 	chi21.val = 0.5,   # Cross immunity, strain 2 against strain 1 (def. 0.3)
# 	amplitude = 0.66,     # R0 seasonality amplitude (dev. 1)
# 	baseline = 1.4,    # R0 seasonality baseline (def. 1.5)
# 	phi.val = -3.8,      # R0 seasonality phase shift, weeks (def. -4)
# 	kappa.val = 0.01,  # Force of infection pulse size for introductions  
# 	importtime1 = 0,  # Time of strain 1 importation, weeks
# 	importtime2 = 52,   # Time of strain 2 importation, weeks
# 	importlength = 0.5 # Duration of importation pulse, weeks
# 	)


# A test with some new values: 
f <- 0.21519566939603102;
maxR0 <- 2.1986663426600286;
parms=c(
	sigma1.val = 1/44.94837857251785, # Waning immunity rate, strain 1, weeks (def. 1/40)
	sigma2.val = 1/44.96099593027429, # Waning immunity rate, strain 2, weeks (def. 1/40)
	nu.val = 1/(2.999751226892689/7),      # Rate of progression to infection, weeks (def. 1/1)
	gamma.val = 1/(5.043512653991342/7),   # Rate of recovery, weeks (def. 1/1)
	chi12.val = 0.7789022405467717,   # Cross immunity, strain 1 against strain 2 (def. 0.7)
	chi21.val = 0.5068431811943369,   # Cross immunity, strain 2 against strain 1 (def. 0.3)
	amplitude = f*maxR0,     # R0 seasonality amplitude (dev. 1)
	baseline = maxR0-f*maxR0,    # R0 seasonality baseline (def. 1.5)
	phi.val = 2.036402516800722,      # R0 seasonality phase shift, weeks (def. -4)
	kappa.val = 0.01,  # Force of infection pulse size for introductions  
	importtime1 = 52,  # Time of strain 1 importation, weeks
	importtime2 = 0,   # Time of strain 2 importation, weeks
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

odefit <- as.data.frame(lsoda(y,times,model,parms, hmax=0.2)) %>% as_tibble()

# From the original paper:
# scalingfactor <- 0.075
# A test:
scalingfactor <- 0.05993353450606306

extractrange <- c(52*23.5, 52*28.5)

figodefit_nrevss <- odefit %>% 
	filter(time>=extractrange[1] & time<=extractrange[2]) %>%
	mutate(Strain1=I1S2 + I1E2 + I1I2 + I1R2) %>%
	mutate(Strain2=S1I2 + E1I2 + I1I2 + R1I2) %>%
	mutate(WEEKINDEX=1:nrow(.)) %>%
	left_join(nrevssCDC_ILI, by="WEEKINDEX") %>%
	select(WEEKEND, Strain1, Strain2, OC43_ILI, HKU1_ILI) %>%
	mutate(OC43_ILI=OC43_ILI) %>%
	mutate(HKU1_ILI=HKU1_ILI) %>% 
	mutate(Strain1=Strain1*scalingfactor) %>% 
	mutate(Strain2=Strain2*scalingfactor) %>% 
	rename("OC43"=OC43_ILI, "HKU1"=HKU1_ILI) %>%
	pivot_longer(c("Strain1","Strain2","OC43","HKU1"), names_to="Virus", values_to="ProportionInfected") %>%
	mutate(DataType=case_when(Virus%in%c("OC43","HKU1")~"Actual",TRUE~"Simulated")) %>%
	mutate(Virus=case_when(Virus=="Strain1"~"OC43", Virus=="Strain2"~"HKU1", TRUE~Virus)) %>%
	ggplot(aes(x=WEEKEND, y=100*ProportionInfected, col=Virus, linetype=DataType)) + 
		geom_line(size=1.2, alpha=0.6) + 
		scale_color_manual(limits=c("OC43","HKU1"), values=c("OC43"="blue", "HKU1"="red")) + 
		labs(x="Year", y="% positive x % ILI") + 
		theme_minimal() + 
		theme(text=element_text(size=20), legend.title=element_blank())
# ggsave(figodefit_nrevss, file="/Users/sk792/DropboxHarvard/Projects/Wuhan/writeup/seasonality/figures/odefit_nrevss.pdf", width=8, height=5)
# ggsave(figodefit_nrevss, file="/Users/sk792/DropboxHarvard/Projects/Wuhan/writeup/seasonality/figures/odefit_nrevss.png", width=8, height=5)

# ggsave(figodefit_nrevss, file="/Users/sk792/DropboxHarvard/Projects/Wuhan/writeup/Science/Submission2/figures/fig2A.pdf", width=8, height=5)
# ggsave(figodefit_nrevss, file="/Users/sk792/DropboxHarvard/Projects/Wuhan/writeup/Science/Submission2/figures/fig2A.png", width=8, height=5)



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


# Add in some opacity: --------------------------------------------------------

alphavals <- nrevssCDC_ILI %>% 
	mutate(OC43_alpha = OC43_ILI/max(.$OC43_ILI)) %>%
	mutate(HKU1_alpha = HKU1_ILI/max(.$HKU1_ILI)) %>%
	select(WEEKEND, OC43_alpha, HKU1_alpha) %>%
	pivot_longer(-WEEKEND, names_to="Virus", values_to="alpha") %>%
	mutate(Virus=case_when(Virus=="OC43_alpha"~"OC43", Virus=="HKU1_alpha"~"HKU1")) 


redf_sim_alpha <- left_join(redf_sim, alphavals, by=c("WEEKEND","Virus"))
redf_wt_alpha <- left_join(redf_wt, alphavals, by=c("WEEKEND","Virus"))

figwtre_alpha <- ggplot() + 
	geom_line(data=redf_sim_alpha, aes(x=WEEKEND, y=Re, col=Virus), size=1.2, alpha=0.6) + 
	geom_point(data=redf_wt_alpha, aes(x=WEEKEND, y=Re, col=Virus, alpha=alpha), size=2) +
	# geom_line(data=redf_wt, aes(x=WEEKEND, y=Re, col=Virus), stat="smooth", method="loess", span=0.2, alpha=0.7, size=1) +
	scale_color_manual(limits=c("OC43","HKU1"), values=c("OC43"="blue", "HKU1"="red")) + 
	labs(x="Week", y="Re") + 
	theme_minimal() + 
	theme(legend.position="none", text=element_text(size=20)) + 
	facet_wrap(~Virus, nrow=1)
# figwtre_alpha
# ggsave(figwtre_alpha, file="/Users/sk792/DropboxHarvard/Projects/Wuhan/writeup/seasonality/figures/wtre_alpha.pdf", width=12, height=5)
# ggsave(figwtre_alpha, file="/Users/sk792/DropboxHarvard/Projects/Wuhan/writeup/seasonality/figures/wtre_alpha.png", width=12, height=5)

# ggsave(figwtre_alpha, file="/Users/sk792/DropboxHarvard/Projects/Wuhan/writeup/Science/Submission2/figures/fig2BC.pdf", width=12, height=5)
# ggsave(figwtre_alpha, file="/Users/sk792/DropboxHarvard/Projects/Wuhan/writeup/Science/Submission2/figures/fig2BC.png", width=12, height=5)


figwtre_size <- ggplot() + 
	geom_line(data=redf_sim_alpha, aes(x=WEEKEND, y=Re, col=Virus), size=1.2, alpha=0.6) + 
	geom_point(data=redf_wt_alpha, aes(x=WEEKEND, y=Re, col=Virus, alpha=alpha, size=alpha)) +
	# geom_line(data=redf_wt, aes(x=WEEKEND, y=Re, col=Virus), stat="smooth", method="loess", span=0.2, alpha=0.7, size=1) +
	scale_color_manual(limits=c("OC43","HKU1"), values=c("OC43"="blue", "HKU1"="orange")) + 
	labs(x="Week", y="Re") + 
	theme_minimal() + 
	theme(legend.position="none", text=element_text(size=20)) + 
	facet_wrap(~Virus, nrow=1)

	
