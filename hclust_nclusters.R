
flag.save = 0
by.conc = 1

by.height = 0
n.cl = 	4# number of clusters


stat_weights = c(1,1,2)


ID.test = 2
cc.test = c(0, 33, 100, 333, 1000)
cc.test = 0

#@ might want to do three clusters with 1000
#6 might want to do three clusters with 333

flag.home = 1

flag.link = 2
linkage.type = c("ward.D2","average", "complete","centroid", "median", "single" )
link.k = linkage.type[flag.link]

library(scatterplot3d)


#mycolors = c("Pink1", "Maroon3",  "orange", "SpringGreen2", "DeepSkyBlue2", "MediumPurple2", "slategray", "Gold")
classcol = c("tomato2", "springgreen3","dodgerblue", "MediumPurple2", "Maroon3", "slategray", "gold") #stuck, sub, free)

mycolors = c( "Maroon3",  "SpringGreen3", "MediumPurple2", "slategray", "Gold")

DonorID = c("F02", "F05", "F08","F13","F15", "F17", "F21")
Concentration = c("c0","c33","c100","c333","c1000")
Conc = c("0","33","100","333","1000")
numeric.conc = c(0, 0.033, 0.1, 0.333, 1.0)

endo_conc = (c(25, 70, 48, 26, 80, 190, 18))/rep(1000, 7) #by donor # UNITS: mu g/mL

TOL = 1
dt = 1/15
pix_conv = 0.229442


hclust.class = c()
means.table.c = list()


per.imm.table = mat.or.vec(length(DonorID), length(Conc))
colnames(per.imm.table) = Conc
rownames(per.imm.table) = DonorID






if( flag.home == 0){
load("/Users/Mjensen1/Dropbox/Ab-Virion-Mucin Modeling and Data/Statistical Analysis/final_data.Rdata")

load("/Users/Mjensen1/Dropbox/Ab-Virion-Mucin Modeling and Data/RData_and_Documentation/probability_weighted_immobilized_all.Rdata")

outliers = c(123,947,1670,3010)-1
mydata = mydata[-outliers,]

}else if(flag.home == 1 ){
	
load("/Users/RachelJensen/Dropbox/Ab-Virion-Mucin Modeling and Data/Statistical Analysis/final_data.Rdata")

load("/Users/RachelJensen/Dropbox/Ab-Virion-Mucin Modeling and Data/RData_and_Documentation/probability_weighted_immobilized_all.Rdata")
#mydata = mydata[-c(1669, 1849 ),]
#outliers = c(123,947,2013,3010)-1
#mydata = mydata[-outliers,]

}


mydata.log = mydata
mydata.log[,7] = log10(mydata[,7])





if(by.conc == 1){
for(ID in ID.test){

donor.data.index = which(mydata[,1] == ID)


fac.data = mydata.log[donor.data.index,]
dim(fac.data)
fac.data$exoAB = as.factor(fac.data$exoAB)

by.exo = split(fac.data, fac.data$exoAB)
lapply(by.exo, summary)


conc.0 = by.exo$"0"
conc.0 = conc.0[,-c(2,3,8,9)]


conc.33 = by.exo$"0.033"
conc.33 = conc.33[,-c(2,3,8,9)]

conc.100 = by.exo$"0.1" 
conc.100 = conc.100[,-c(2,3,8,9)]

conc.333= by.exo$"0.333"
conc.333 = conc.333[,-c(2,3,8,9)]

conc.1000 = by.exo$"1"
conc.1000 = conc.1000[,-c(2,3,8,9)]

##Hierachical clustering k =  by concentration
scaled.0 = scale(conc.0[,-c(1,2,6)])
w.scaled.0 = scaled.0*matrix(rep(sqrt(stat_weights), times = dim(scaled.0)[1]), ncol = 3, byrow = TRUE)

scaled.33 = scale(conc.33[,-c(1,2,6)])
w.scaled.33 = scaled.33*matrix(rep(sqrt(stat_weights), times = dim(scaled.33)[1]), ncol = 3, byrow = TRUE)


scaled.100 = scale(conc.100[,-c(1,2,6)])
w.scaled.100 = scaled.100*matrix(rep(sqrt(stat_weights), times = dim(scaled.100)[1]), ncol = 3, byrow = TRUE)

scaled.333 = scale(conc.333[,-c(1,2,6)])
w.scaled.333 = scaled.333*matrix(rep(sqrt(stat_weights), times = dim(scaled.333)[1]), ncol = 3, byrow = TRUE)

scaled.1000 = scale(conc.1000[,-c(1,2,6)])
w.scaled.1000 = scaled.1000*matrix(rep(sqrt(stat_weights), times = dim(scaled.1000)[1]), ncol = 3, byrow = TRUE)


par(mfcol= c(2,5), mar = c(4,4,2,2), oma = c(1,1,1,1))
c.index = 0
for(cc.loop in cc.test){
c.index = c.index + 1
flag.c = cc.loop

 if(ID == 1 && c.index == 5){
	n.cl = 5
}else if(ID ==4 ){n.cl = 5
	}else{n.cl = 4}
#-------------------------------
# Exogenous = 0
#-------------------------------
if (flag.c == 0){

dis.matrix.0 = dist(w.scaled.0, method = "euclidean")

h.ward.0 = hclust(dis.matrix.0, method = link.k)

#complete, centroid, median

h.ward = h.ward.0
c.data = conc.0

}else if(flag.c == 33){
#-------------------------------
# Exogenous = 33
#-------------------------------

dis.matrix.33 = dist(w.scaled.33, method = "euclidean")
h.ward.33 = hclust(dis.matrix.33, method = link.k)

h.ward = h.ward.33
c.data = conc.33

#plot(h.ward.33)


}else if(flag.c == 100){
#-------------------------------
# Exogenous = 100
#-------------------------------

dis.matrix.100 = dist(w.scaled.100, method = "euclidean")
h.ward.100 = hclust(dis.matrix.100, method = link.k)
#plot(h.ward.100)
h.ward = h.ward.100
c.data = conc.100


}else if(flag.c == 333){
#-------------------------------
# Exogenous = 333
#-------------------------------

dis.matrix.333 = dist(w.scaled.333, method = "euclidean")

h.ward.333 = hclust(dis.matrix.333, method = link.k)
h.ward = h.ward.333
c.data = conc.333


}else if(flag.c == 1000){
#-------------------------------
# Exogenous = 1000
#-------------------------------

dis.matrix.1000 = dist(w.scaled.1000, method = "euclidean")

h.ward.1000 = hclust(dis.matrix.1000, method = link.k)
h.ward = h.ward.1000
c.data = conc.1000
}else{
	print("pick from (0,33, 100, 333, 1000)")
}


#-------------------------------
# Constructing the tree
#-------------------------------

plot(h.ward, main = paste("Concentration", flag.c), xlab = " ")


if(by.height == 0){
rect.ward = rect.hclust(h.ward, k = n.cl)
}else{
rect.ward = rect.hclust(h.ward, h = set_height)
n.cl = length(rect.ward)
}


group_k = list()
for(k in 1:n.cl){
	group_k[[k]] = unique(rect.ward[[k]])
}

gk_data = list()

for(k in 1:n.cl){
	gk_data[[k]] = c.data[group_k[[k]],]
}




means.table.c[[c.index]] = mat.or.vec(n.cl ,3)

for(k in 1:n.cl){
means.table.c[[c.index]][k,] = apply(gk_data[[k]][,-c(1,2,6)], 2, mean)
}
colnames(means.table.c[[c.index]]) = c("xACF", "yACF", " Diffusivity")
rownames(means.table.c[[c.index]]) = paste("Hclust", 1:n.cl)


order.diff = order(means.table.c[[c.index]][,3])

group_k = group_k[order.diff]
gk_data = gk_data[order.diff]

cluster.class = mat.or.vec(nr = dim(c.data)[1], nc = 1)

for(k in 1:n.cl){
cluster.class[group_k[[k]]] = k 

}	
hclust.class = c(hclust.class, cluster.class)




diff.label = c(expression(10^-3), expression(10^-2), expression(10^-1),expression(10^0),expression(10^1))

 plot(gk_data[[1]][,5], gk_data[[1]][,3], col = classcol[1], pch =4, cex = 0.8, ylim = c(-1,1), ylab = expression(paste("ACF (h =1)")), xlab =expression(paste("Estimated Diffusivity (", mu, m^2, s^-1, ")" )), yaxt = "n", mar = c(0.5,0.5,0.5,0.5), xaxt = "n", xlim = c(-3.2,1))
 points(gk_data[[1]][,5], gk_data[[1]][,4], col = classcol[1]  ,pch = 6, cex = 0.8)

 
 axis(2, at = seq(-1,1,by = 0.5), label =seq(-1,1,by = 0.5), las = 2  )
 axis(1, at = seq(-3,1,by =1), label =diff.label )
	
for(k in 2:n.cl){
points(gk_data[[k]][,5], gk_data[[k]][,3], col = classcol[k]  ,pch = 4, cex = 0.8)
points(gk_data[[k]][,5], gk_data[[k]][,4], col = classcol[k]  ,pch = 6, cex = 0.8)

}
#donuts = c(1,2,4, 5,9,13, 16,17, 19, 21)
#points(gk_data[[k]][donuts,5], gk_data[[k]][donuts,3] ,pch = 4, cex = 0.8)


	min_CI = -1*qnorm((1+0.95)/2)/sqrt(max(c.data[,2]))
	max_CI = -1*qnorm((1+0.95)/2)/sqrt(min(c.data[,2]))


	abline(0,0)
	abline(-max_CI, 0, lty = 3, col = "slategray");abline(max_CI, 0, lty = 3, col = "slategray")
	abline(-min_CI, 0, lty = 3, col = "slategray");abline(min_CI, 0, lty = 3, col = "slategray")






} ##End of concentration loop

title(paste(" Donor ",DonorID[ID] ) , outer = TRUE)

} #End of the ID loop
} #Clustering by concentration











if(by.conc == 0){

for(ID in ID.test){

donor.data.index = which(mydata[,1] == ID)


donor.data = mydata.log[donor.data.index, ]
scaled.data = scale(donor.data[,5:7])


w.scaled.data = scaled.data*matrix(rep(sqrt(stat_weights), times = dim(donor.data)[1]), ncol = 3, byrow = TRUE)



dis.matrix = dist(w.scaled.data, method = "euclidean")

h.ward = hclust(dis.matrix, method = link.k)


#-------------------------------
# Constructing the tree
#-------------------------------
par(mfrow= c(2,3), mar = c(4,4,2,2), oma = c(1,1,1,1))

plot(h.ward, main = paste("Dendogram"), xlab = " ")


if(by.height == 0){
rect.ward = rect.hclust(h.ward, k = n.cl)
}else{
rect.ward = rect.hclust(h.ward, h = set_height)
n.cl = length(rect.ward)
}


group_k = list()
for(k in 1:n.cl){
	group_k[[k]] = unique(rect.ward[[k]])
}

gk_data = list()

for(k in 1:n.cl){
	gk_data[[k]] = donor.data[group_k[[k]],]
}




means.table = mat.or.vec(n.cl ,3)

for(k in 1:n.cl){
means.table[k,] = apply(gk_data[[k]][,5:7], 2, mean)
}
colnames(means.table) = c("xACF", "yACF", " Diffusivity")
rownames(means.table) = paste("Hclust", 1:n.cl)

if(means.table[2,3] > means.table[1,3]){
	hclust.class[donor.data.index[group_k[[2]]]]  = 1
	classcol = c("tomato2", "aquamarine3","dodgerblue", "MediumPurple2", "Maroon3", "slategray", "gold") 
	
}else{
	hclust.class[donor.data.index[group_k[[1]]]]  = 1		
	classcol = c("aquamarine3", "tomato2","dodgerblue", "MediumPurple2", "Maroon3", "slategray", "gold") 
}




diff.label = c(expression(10^-3), expression(10^-2), expression(10^-1),expression(10^0),expression(10^1))





for(cc in 1:5){

	conc.index = which(donor.data$exoAB == numeric.conc[cc])	
	conc.group.index = lapply(group_k, intersect, conc.index)
	
	cg.1 = c(conc.group.index[[1]])

	 plot(donor.data[cg.1,7], donor.data[cg.1,5], col = classcol[1], pch =4, cex = 0.8, ylim = c(-1,1), ylab = expression(paste("ACF (h =1)")), xlab =expression(paste("Estimated Diffusivity (", mu, m^2, s^-1, ")" )), yaxt = "n", mar = c(0.5,0.5,0.5,0.5), xaxt = "n", xlim = c(-3.2,1))
 	points(donor.data[cg.1,7], donor.data[cg.1,6], col = classcol[1], pch =6, cex = 0.8, ylim = c(-1,1), ylab = expression(paste("ACF (h =1)")), xlab =expression(paste("Estimated Diffusivity (", mu, m^2, s^-1, ")" )), yaxt = "n", mar = c(0.5,0.5,0.5,0.5), xaxt = "n", xlim = c(-3.2,1))

 
 
 
 
	 axis(2, at = seq(-1,1,by = 0.5), label =seq(-1,1,by = 0.5), las = 2  )
	 axis(1, at = seq(-3,1,by =1), label =diff.label )
	title(paste("Concentration", numeric.conc[cc]))
	for(k in 2:n.cl){
		cg.k = c(conc.group.index[[k]])
		points(donor.data[cg.k,7], donor.data[cg.k,5], col = classcol[k]  ,pch = 4, cex = 0.8)
		points(donor.data[cg.k,7], donor.data[cg.k,6], col = classcol[k]  ,pch = 6, cex = 0.8)

	}


	

	min_CI = -1*qnorm((1+0.95)/2)/sqrt(max(donor.data[,4]))
	max_CI = -1*qnorm((1+0.95)/2)/sqrt(min(donor.data[,4]))


	abline(0,0)
	abline(-max_CI, 0, lty = 3, col = "navy");abline(max_CI, 0, lty = 3, col = "navy")
	abline(-min_CI, 0, lty = 3, col = "navy");abline(min_CI, 0, lty = 3, col = "navy")






} ##End of concentration loop

title(paste(" Donor ",DonorID[ID] ) , outer = TRUE)




} #End of the ID loop


} # Not clustering by concentration









#if(flag.save == 1){
	#data.hclust.class = cbind(mydata, hclust.class)
	#hclust.percent.imm =  per.imm.table 
	#save(data.hclust.class, hclust.percent.imm, file ="data_hclust_final_prune.Rdata")
	
#}

	#data.hclust.class = cbind(mydata, hclust.class)
	#hclust.percent.imm =  per.imm.table 
	#save(data.hclust.class,  file ="data_hclust_final_4class.Rdata")
	




## For the figure
 #class.fig.data = cbind(mydata[donor.data.index,], hclust.class)
#save(class.fig.data, file = "F08_data_4classification.Rdata")

