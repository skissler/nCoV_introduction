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
# Internal functions
# =============================================================================

beta.val <- function(t, amplitude=1, baseline=1.5, phi.val=-4, gamma.val=1){
	gamma.val*(amplitude/2 * cos(2*pi*(t-phi.val)/52) + (amplitude/2 + baseline))
}

p1 <- function(t,kappa.val,importtime1,importlength){
	ifelse(t>importtime1 & t<=(importtime1+importlength), kappa.val, 0)
}

p2 <- function(t,kappa.val,importtime2,importlength){
	ifelse(t>importtime2 & t<=(importtime2+importlength), kappa.val, 0)
}

# =============================================================================
# State update functions 
# =============================================================================

S1S2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	beta.val(t,amplitude,baseline,phi.val,gamma.val)*(I1S2+I1E2+I1I2+I1R2)*S1S2 + p1(t,kappa.val,importtime1,importlength)*S1S2
	})}
S1S2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	beta.val(t,amplitude,baseline,phi.val,gamma.val)*(S1I2+E1I2+I1I2+R1I2)*S1S2 + p2(t,kappa.val,importtime2,importlength)*S1S2
	})}
E1S2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	nu.val*E1S2
	})}
E1S2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	(1-chi12.val)*beta.val(t,amplitude,baseline,phi.val,gamma.val)*(S1I2+E1I2+I1I2+R1I2)*E1S2 + (1-chi12.val)*p2(t,kappa.val,importtime2,importlength)*E1S2
	})}
S1E2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	(1-chi21.val)*beta.val(t,amplitude,baseline,phi.val,gamma.val)*(I1S2+I1E2+I1I2+I1R2)*S1E2 + (1-chi21.val)*p1(t,kappa.val,importtime1,importlength)*S1E2
	})}
S1E2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	nu.val*S1E2
	})}
E1E2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	nu.val*E1E2 
	})}
E1E2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	nu.val*E1E2 
	})}
I1S2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	gamma.val*I1S2
	})}
I1S2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	(1-chi12.val)*beta.val(t,amplitude,baseline,phi.val,gamma.val)*(S1I2+E1I2+I1I2+R1I2)*I1S2 + (1-chi12.val)*p2(t,kappa.val,importtime2,importlength)*I1S2
	})}
S1I2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	(1-chi21.val)*beta.val(t,amplitude,baseline,phi.val,gamma.val)*(I1S2+I1E2+I1I2+I1R2)*S1I2 + (1-chi21.val)*p1(t,kappa.val,importtime1,importlength)*S1I2
	})}
S1I2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	gamma.val*S1I2
	})}
R1S2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	sigma1.val*R1S2
	})}
R1S2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	(1-chi12.val)*beta.val(t,amplitude,baseline,phi.val,gamma.val)*(S1I2+E1I2+I1I2+R1I2)*R1S2 + (1-chi12.val)*p2(t,kappa.val,importtime2,importlength)*R1S2
	})}
I1E2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	gamma.val*I1E2
	})}
I1E2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	nu.val*I1E2
	})}
E1I2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	nu.val*E1I2
	})}
E1I2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	gamma.val*E1I2
	})}
S1R2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	(1-chi21.val)*beta.val(t,amplitude,baseline,phi.val,gamma.val)*(I1S2+I1E2+I1I2+I1R2)*S1R2 + (1-chi21.val)*p1(t,kappa.val,importtime1,importlength)*S1R2
	})}
S1R2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	sigma2.val*S1R2
	})}
R1E2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	sigma1.val*R1E2
	})}
R1E2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	nu.val*R1E2
	})}
I1I2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	gamma.val*I1I2 
	})}
I1I2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	gamma.val*I1I2
	})}
E1R2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	nu.val*E1R2 
	})}
E1R2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	sigma2.val*E1R2
	})}
R1I2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	sigma1.val*R1I2
	})}
R1I2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	gamma.val*R1I2
	})}
I1R2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	gamma.val*I1R2
	})}
I1R2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	sigma2.val*I1R2
	})}
R1R2c1 <- function(t,y,parms){with(as.list(c(y,parms)),{
	sigma1.val*R1R2
	})}
R1R2c2 <- function(t,y,parms){with(as.list(c(y,parms)),{
	sigma2.val*R1R2
	})}
