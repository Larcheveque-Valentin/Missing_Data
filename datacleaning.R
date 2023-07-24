
setwd("C:/Users/Utilisateur/Documents/Master SSD/Donnéess manquantes/Data-cleaning")

prenom_val=read.table("PRENOMS VALIDES.csv",header=T ,stringsAsFactors = F, sep="\t", dec="." )
prenom_val=as.data.frame( gsub(pattern = "\x83",x=prenom_val[,1],replacement = "E"))
prenom_val=as.data.frame( gsub(pattern = "\xe8",x=prenom_val[,1],replacement = "Ê"))
colnames(prenom_val)<-"PRENOMS VALIDES"
prenom_val

prenom_a_verif=read.table("PRENOMS A VERIFIER.csv",header=T ,stringsAsFactors = F, sep="\t", dec="." )
prenom_a_verif=as.data.frame( gsub(pattern = "\x83",x=prenom_a_verif[,1],replacement = "E"))
prenom_a_verif=as.data.frame( gsub(pattern = "\xe8",x=prenom_a_verif[,1],replacement = "Ê"))
colnames(prenom_a_verif)<-"PRENOMS A VERIFIER"
prenom_a_verif



library("stringdist")
Corrector<- function(true_names, to_be_corrected){
  false_names=to_be_corrected
      stringd<-matrix(0,dim(true_names)[1],dim(to_be_corrected)[1])
for ( i in 1:dim(true_names)[1]){
  for (j in 1:dim(to_be_corrected)[1]){
      stringd[i,j]<-stringdist(true_names[i,],to_be_corrected[j,],"lv")
  }
}
      
for (j in 1:dim(to_be_corrected)[1]){
  if(min(stringd[,j])!=0){
        to_be_corrected[j,]<-true_names[which.min(stringd[,j]),]
        
  }}
      diff<-to_be_corrected!=false_names
      corrected_names<-data.frame(to_be_corrected,diff)
      colnames(corrected_names)<- c("Corrected names", "Was corrected:")
    return(corrected_names) 
}
system.time(prenom_corrected<-Corrector(true_names = prenom_val,to_be_corrected = prenom_a_verif))
prenom_corrected


