#' Seeded fast causal inference algorithm
#'
#' Uses information from the domain problem encoded in known paths of a partial structure
#' 
#' @param G.am An adjacency matrix of the initial CPDAG containig one or more undefined causal edges to be resolved. The edgemark-code in am refers to the row index.
#' @param V.data A matrix with \code{"m"} observations (rows) of the \code{"p"} variables (columns)
#'
#' @return A adjacency matrix with the found equivalence class represented by a partial ancestral graph
#' encoded in a csv file.
#' @author Samuel Montero-Hernandez, \email{samuel@@inaoep.mx}
#' @keywords graphical models

rm(list = ls())
library("pcalg")
library("Rgraphviz")
source("skeletonvCondS.r")

# Example of use with predefined connectome links retrieved from Hagmann et al. 2008 PLoS Biol 6(7):e159. and Joshi et al. 2010. IEEE International Symposium on Biomedical Imaging 844–847.

chlabels = as.character(1:24)
elem <- 24*24
fixEd <-  matrix(rep(FALSE,elem),nr=24,nc=24)
fixEd[14,10]<-TRUE
fixEd[14,13]<-TRUE
fixEd[14,2]<-TRUE
fixEd[14,1]<-TRUE
fixEd[14,20]<-TRUE

fixEd[1,2]<-TRUE
fixEd[1,13]<-TRUE
fixEd[1,14]<-TRUE
fixEd[1,20]<-TRUE

fixEd[20,10]<-TRUE
fixEd[20,2]<-TRUE
fixEd[20,1]<-TRUE
fixEd[20,14]<-TRUE
fixEd[20,13]<-TRUE

fixEd[10,20]<-TRUE
fixEd[10,13]<-TRUE
fixEd[10,14]<-TRUE
fixEd[10,2]<-TRUE

fixEd[13,20]<-TRUE
fixEd[13,10]<-TRUE
fixEd[13,14]<-TRUE
fixEd[13,1]<-TRUE
fixEd[13,2]<-TRUE

fixEd[2,10]<-TRUE
fixEd[2,20]<-TRUE
fixEd[2,1]<-TRUE
fixEd[2,14]<-TRUE
fixEd[2,13]<-TRUE
#load fnirs data
#"vts_HbO2_nov.csv"
#"vts_HbO2_train.csv"
#"vts_HbO2_exp.csv"
#"vts_HbO2_all.csv"
pags<-c()
group<-c("nov","train","exp","all")
i<-1
fixGap <- !fixEd
for(gp in group){
	channels <- read.csv(file=paste("vts_HbO2_",gp,".csv",sep=""),head=TRUE,sep=",")
	suffStat <- list(C = cor(channels), n = nrow(channels))
	#tmp <- rfci(suffStat, indepTest = gaussCItest, alpha = 0.05, labels = chlabels)
	#pags<-c(pags,tmp)
	tmp <- rfci(suffStat, indepTest = gaussCItest, alpha = 0.05, labels = chlabels, fixedGaps=fixGap)
	pags<-c(pags,tmp)

	 par(mfrow = c(1, 1))
	 plot(pags[[i]])
	 #title(sub="No Connectome")
	 #plot(pags[[i+1]])
	 title(sub="With Connectome")


	write.table(pags[[i]]@amat, paste("adjmatFixGaps",gp,"C.csv",sep=""),sep=",")
	#write.table(pags[[i+1]]@amat, paste("adjmatGaps",gp,"C.csv",sep=""),sep=",")
	i<-i+1
}