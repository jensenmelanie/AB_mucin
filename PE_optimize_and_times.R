

## File Summary
#2D parameter space data

at.home = 0
flag.save = 1
##---------------------------------------------------------------------
## Parameters to pick
##---------------------------------------------------------------------
Thresh.vec = 1:60
##---------------------------------------------------------------------
## Assumed parameter values
##---------------------------------------------------------------------


knockdown_factor = 0.9
moff = 100
kd = 0.8969
monM = (1 - knockdown_factor)/knockdown_factor * moff
#monM = ((1- knockdown_factor)/(knockdown_factor))*0.01
Ab.tested = c(0, 0.0333, 0.1, 0.333, 1)

#_________________________________________________________________________

## Estimated Paramters


Nstar = c(65,seq(70,495,by = 5), seq(500,1000, by = 10))
#mod.c.log = seq(0,2 ,length = 100)
#mod.c = c(10^mod.c.log, 200, 300, 500)

mod.c = unique(10^c(seq(-2, 2, by = 0.05),seq(2, 4, by = 0.025) ))


p_int = seq(0.5, 0.999, by = 0.005)#length = 50)
endo = c(0.0001,seq(0.001, 0.01, by = 0.001), seq(0.01, 0.25, by= 0.0025))


# # 
   ## Nstar = seq(100, 1000, by = 100)
   ## mod.c.log = seq(0, 4, length = 10)
   ## mod.c = 10^mod.c.log
   ## p_int = seq(0.5, 0.99, length = 10)
   ##endo = sort(c(0.0291,seq(0.001, 0.25, length = 10))) #20


#______________________________________________________________________________

## Percent Immobilized Tables 
#______________________________________________________________________________



Conc = c("0",".033",".10",".333","1.00")
DonorID = c("F02", "F05", "F08","F13","F15", "F17", "F21")


if( at.home == 1){
	setwd("/Users/RachelJensen/Dropbox/Ab-Virion-Mucin Modeling and Data/Final_model_fitting")
	load("/Users/RachelJensen/Dropbox/Ab-Virion-Mucin Modeling and Data/November_2017_codes/data_hclust_final_prune2.Rdata")
	
	load("/Users/RachelJensen/Dropbox/Ab-Virion-Mucin Modeling and Data/Final_model_fitting/final_p_imm.Rdata")

		class_stuck_all = weighted.p.imm

}else{
		setwd("/Users/Mjensen1/Dropbox/Ab-Virion-Mucin Modeling and Data/Final_model_fitting")

		load("/Users/Mjensen1/Dropbox/Ab-Virion-Mucin Modeling and Data/November_2017_codes/data_hclust_final_prune2.Rdata")

	
		load("/Users/Mjensen1/Dropbox/Ab-Virion-Mucin Modeling and Data/Final_model_fitting/final_p_imm.Rdata")

	
	class_stuck_all = weighted.p.imm
}



num.particles = matrix(nrow = 7, ncol = 5) 


for(z in 1:7){
	for(cc in 1:5){
		num.particles[z,cc] = length(intersect(which(data.hclust.class$DonorID == z), which(data.hclust.class$exoAB == unique(data.hclust.class$exoAB)[cc])))
	}
	
}




##----------------------------------
##Required function
##----------------------------------

E.stuck.fnc = function(Nstar, mod.c, exo, endo, Thresh){
	
	total.ab = endo + exo

n.vec = Thresh:Nstar

p.av = total.ab/(total.ab + kd)
p.mav.c = mod.c*monM/(moff + mod.c*monM)


summands = (1/(Thresh*moff))*(1- unlist(lapply(as.list(n.vec), pbinom, q = (Thresh -1), prob = p.mav.c)))*(choose(Nstar, n.vec)/choose(n.vec, Thresh))*(kd/(kd+total.ab))^(Nstar-n.vec)*(1/(p.mav.c)^Thresh)*(((total.ab*(moff+ mod.c*monM))/((total.ab + kd)*moff))^n.vec)*(moff/(moff+ mod.c*monM))^Thresh


E.stuck.fnc = sum(summands)
	
	

}



E.free.fnc <- function(Nstar,mod.c, exo, endo, Thresh){

	total.ab = endo + exo

	n.vec = Thresh:Nstar
	
	pi.Tminus = unlist(lapply(as.list(n.vec), dbinom, x = (Thresh-1), prob = monM/(monM+moff)))
	
	prob.n = dbinom(n.vec, Nstar, total.ab/(total.ab + kd))
	
	E.n = (1/((n.vec -Thresh +1)*monM))*(unlist(lapply(as.list(n.vec), pbinom, q = (Thresh-1), prob = monM/(monM + moff))))/pi.Tminus

	

 	E.free.fnc = sum(E.n*prob.n)	


}



##----------------------------------
## Computing the Expected time functions
##----------------------------------
tested.exo = c(0, 0.0333, 0.1, 0.333, 1)

frac.imm.n.Thresh = list()

for(tt in 1:length(Thresh.vec)){
	
Thresh = Thresh.vec[tt]	
	

n.vec = Thresh:max(Nstar)
num.n = max(Nstar)-Thresh + 1


pi.Tminus = unlist(lapply(as.list(n.vec), dbinom, x = (Thresh-1), prob = monM/(monM+moff)))


E.tau.vec = (1/((n.vec -Thresh +1)*monM))*(unlist(lapply(as.list(n.vec), pbinom, q = (Thresh-1), prob = monM/(monM + moff))))/pi.Tminus



## -------Get the E.sigma.vector for all n from Thresh:Nstar----
## Rows are with differenct modc values, columns are different n

E.sigma.mat = mat.or.vec(nr = length(mod.c), nc = num.n)
frac.imm.n.Thresh[[tt]] = mat.or.vec(nr = length(mod.c), nc = num.n)

for( cc in 1:length(mod.c)){
p.mav.c = (mod.c[cc]*monM)/(mod.c[cc]*monM + moff)


pi.T = unlist(lapply(as.list(n.vec), dbinom, x = Thresh, prob = p.mav.c))


E.sigma.mat[cc,] = (1/(Thresh*moff))*(1- unlist(lapply(as.list(n.vec), pbinom, q = (Thresh -1), prob = p.mav.c)))/pi.T


E.sigma.mat[cc,which(E.sigma.mat[cc,] == "Inf")] = max(E.sigma.mat[1:(cc-1),])


frac.imm.n.Thresh[[tt]][cc,] = E.sigma.mat[cc,]/(E.sigma.mat[cc,] + E.tau.vec)


}

}




##----------------------------------
## Start of the donor by Donor loop
##----------------------------------

test.ID = 6:7

for(ID in 1:length(DonorID)){
#for(ID in test.ID){

##---------------------------------------------------------------------
## Output variables
##---------------------------------------------------------------------
min.observed.table = mat.or.vec(length(Thresh.vec), 5)
colnames(min.observed.table) = c("c", "Nstar", "q", "endo", "chi")

residual.matrices = list()
opt.endo.matrices = list()
opt.q.matrices = list()



observed = mat.or.vec(length(ID),5)
observed[1:length(ID),] = class_stuck_all[ID,]
num.particle.ID = num.particles[ID,]

DonorID = c("F02", "F05", "F08","F13","F15", "F17", "F21")
line_ID = c(1,2,1,4,5,6,2)
points_ID = c(0, 1, 17, 4, 8, 15, 6)
colors_ID = c("Pink1", "Maroon3", "orange", "SpringGreen2", "DeepSkyBlue2", "MediumPurple2", "slategray")


DonorID = DonorID[ID]
points_ID = points_ID[ID]
line_ID = line_ID[ID]
Color_ID = colors_ID[ID]




##=========================================
##Start of Optimization
##=========================================


ptm <-proc.time()





## -------Stationary dist of Antibody-mucin ----


for(tt in 1:length(Thresh.vec)){
Thresh = Thresh.vec[tt]


prob.imm.almost = function(prob.n, c.value, nstar){
	sumto = nstar- Thresh +1
	sum(frac.imm.n.Thresh[[tt]][c.value,(1:sumto)]*prob.n)	
}



get.res = function(p.interact, prob.imm.al){
	model = p.interact*prob.imm.al
	summands = num.particle.ID*(observed -model)^2/(model*(1-model))
	get.res = sum(summands)
}



		
Z_cNstar = mat.or.vec(nr = length(mod.c), nc = length(Nstar))
overall.opt.endo.index = mat.or.vec(nr = length(mod.c), nc = length(Nstar))
overall.opt.q.index = mat.or.vec(nr = length(mod.c), nc = length(Nstar))

for(cc in 1:length(mod.c)){
	for(ns in 1:length(Nstar)){
	
		min.endo = c()
		opt.q.index = c()

		for(aa in 1:length(endo)){
			total.ab = (tested.exo + endo[aa])
			p.av = total.ab/(total.ab + kd)
			p_N = lapply(as.list(p.av), dbinom, x = (Thresh:Nstar[ns]), size = Nstar[ns])  
			#Each list is for a different concetnration
			
			prob.imm.almost.vec = unlist(lapply(p_N, prob.imm.almost, c.value = cc, nstar = Nstar[ns]))
			prob.imm.almost.mat = matrix(prob.imm.almost.vec, nrow = dim(observed)[1], ncol = 5, byrow = TRUE)	
									
			residuals = unlist(lapply(as.list(p_int),get.res,prob.imm.al = prob.imm.almost.mat ))
		
		 min.endo = c(min.endo, min(residuals))	 		
	 	opt.q.index = c(opt.q.index, which.min(residuals))	
					
		} # end of endo
		
	
	Z_cNstar[cc,ns] = min(min.endo) 
	overall.opt.endo.index[cc,ns] = max(which.min(min.endo),0)
	overall.opt.q.index[cc,ns] = max(opt.q.index[which.min(min.endo)], 0)	
		
	}# end of Nstar

} #End of c/ optimization optimization



Z_tilde = Z_cNstar
Z_tilde[which(Z_cNstar > 7777777777 )] = 7777777777


residual.matrices[[tt]] = Z_tilde
opt.endo.matrices[[tt]] = matrix(endo[overall.opt.endo.index], nrow = dim(overall.opt.endo.index)[1], ncol = dim(overall.opt.endo.index)[2])
opt.q.matrices[[tt]] = matrix(p_int[overall.opt.q.index], nrow = dim(overall.opt.q.index)[1], ncol = dim(overall.opt.q.index)[2])

##----------------------------------
## Levels for the residual 
##----------------------------------

level_log = c((round(log10(min(Z_tilde)),1) -0.1) + 0:3,log10(7777777777)) 


log_Z = matrix(cut(log10(Z_tilde), level_log, include.lowest= TRUE, labels = FALSE), nrow = length(mod.c), ncol = length(Nstar))


chi_gray = gray(seq(1:length(level_log))/length(level_log)) #Just the grayscale 
labelZ = level_log
 
 xxlabel = c(expression(10^-3), expression(10^-2), expression(10^-1),expression(10^0),expression(10^1), expression(10^2))

 xx = log10(mod.c); xx_label = "Mod c"
 yy = Nstar; yy_label = expression(paste(N["*"]) )
 zz = log_Z

 filled.contour(xx, yy, zz, levels = 1:length(level_log), col = chi_gray, ylab = yy_label, xlab = xx_label, ylim = range(yy), xlim= range(xx),  plot.axes ={contour(xx,yy,zz, levels = seq(1,length(level_log)), drawlabels = FALSE, axes = FALSE, frame.plot = FALSE, add = TRUE, col = "black"); axis(1, at = -3:2, label = xxlabel);axis(2)},key.title = title(main = expression(paste(chi[PL]^2))) ,  key.axes = {axis(4, at = seq(1,length(level_log)), labels = labelZ, cex = 0.5)})
 #title(expression(paste("mod.c vs ",N["*"]," vs  ", chi[PL]^2 )))
title(paste("T = ", Thresh))


}# End of Thresh loop





if(flag.save == 1){
filename = paste("PE_", DonorID,".Rdata" , sep = "")
save(residual.matrices, opt.endo.matrices, opt.q.matrices, Nstar, mod.c, p_int, endo,Thresh.vec, file = filename )
}

print(proc.time()- ptm) 
print(DonorID)
} #End Donor loop
