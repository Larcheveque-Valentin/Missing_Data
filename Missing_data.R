
library(norm)
library(missMethods)
library("missForest")


missForest()


library("zoo")
library(dplyr)
library("Rcpp")
library("RcppParallel")
setwd("C:/Users/Utilisateur/Documents/Master SSD/Donnéess manquantes")
porto<-read.table("porto.csv",header = T,sep="\t",stringsAsFactors = T)
porto<-na.omit(porto)



tableau_MAR<-function(dtf,var_exp,var_MAR,dep){
  df<-cbind(matrix(1,length(dtf[,1]),1),dtf[,var_exp])
  qtt<-as.matrix(df)%*%as.matrix(dep)
  qtt<-qtt/sqrt(sum(qtt^2))##calcul de la variable d'int?r?t
  n<-length(df[,1])
  M<-rep(0,n)
  vec_pro<-rep(0,n)
    vec_pro<-exp(qtt)/(1+exp(qtt))
    ## calcul du quantile emprique de
    ##la variable d'int?r?t de l'individu i
    M<-rbinom(n=length(vec_pro),size=1,prob = vec_pro)##tirage au sort du manquemant de x(i,j)
    ##par une loi binomiale param?tr?e par le quantile de la variable d'int?r?t.
    dtf[M==1,var_MAR]<-NA
 
  return(dtf)
}



  tableau_MNAR<-function(df,var_MNAR,dep=c(1,1)){
    dtf<-tableau_MAR(df,var_MNAR,var_MNAR,dep)
    return(dtf)
  }
    
    MCAR<-function(df,p=0.1){
      les_malchanceux_disparues<-rbinom(dim(df)[1]*dim(df)[2],prob = p,size = 1 )
      les_malchanceux_disparues<-matrix(data = les_malchanceux_disparues,
                                        nrow = dim(df)[1],
                                        ncol = dim(df)[2])
      df[les_malchanceux_disparues==1]<-NA
      return(df)}
      
      
    imput_moy <- function(df){
      for (j in 1:ncol(df)){
        vect <- which(is.na(df[j]), arr.ind = TRUE)[, 1]
        if(length(vect) != 0){
          moy <- mean(df[, j], na.rm = TRUE)
          df[vect, j] <- moy
        }
      }
      return(df)
    }
 
    
  
    
    impute_reg <- function(df, col_na, col_reg){
      # col_na sont les colonnes avec des NA
      # col_reg sont les colonnes sur lesquelles on r?gresse, PAS DE NA
      for(j in col_na){
        # Le vecteur des indices des donn?es manquantes
        vect <- which(is.na(df[j]), arr.ind = TRUE)[,1]
        # Les vecteurs des r?gresseurs X et de la variable ? expliquer Y aux indices o? il n'y a pas de NA dans Y
        df_temp <- as.data.frame(cbind(df[-vect, j], df[-vect, col_reg]))
        # On r?gresse Y sur X
        model <- lm(as.matrix(df_temp[1]) ~ as.matrix(df_temp[,2:ncol(df_temp)]), data = df_temp)
        # On ajoute l'intercept
        new_data <- cbind(rep(1, length(vect)), df[vect, col_reg])
        new_data <- as.matrix(new_data)
        # On pr?dit les NA de Y avec le mod?le supra
        # Predict ne marche pas, probl?mes avec le nom des colonnes, donc \hat{Y} = X*\hat{\beta}
        df[vect, j] <- new_data%*%(as.vector(model$coefficients))
      }
      return(df)
      
    }
   
    missSVD<-function(t,rg_approx,n_iter){
     
       NA_places<-is.na(t)
      t<-imput_moy(t)
      t<-as.matrix(t)
      for (i in 1:n_iter){
     X<-t(t)%*%t
    
     Dlambda<-diag(eigen(X)$values)
     U<-eigen(X)$vectors
     W<-(t%*%U)%*%solve(sqrt(Dlambda))
        df<-t(U[,1:rg_approx]%*%sqrt(diag(eigen(X)$values[1:rg_approx]))%*%t(W[,1:rg_approx]))
        
        t<-(t*!NA_places)+(df*NA_places)
      }
      return(t)
    }
    
    
    setClass(Class="mcdudu",slots = c(list_mse="numeric",list_nna="numeric",msetot="numeric",PRESS="numeric"))
  
    Impute_perf_MC<-function(dtf,k=10,method="SVD",maxit_method=100,rg_SVD=2,miss_data="MCAR",p=0.3,ntree=20){
    dtf_test<-dtf  
    list_MSE<-rep(0,k)
      list_nbna<-rep(0,k)
      
      l<-length(dtf[1,])
      q<-length(dtf[,1])
      vect<-rep(1,(l+1))
      if(miss_data=="MCAR"){
        if (method=="miss_forest"){
          for (i in 1:k){
  
            val_manq=rep(0,l)
            
            while(length(val_manq)==0|length(val_manq)==l){
              
              tab<-MCAR(dtf,p)
              val_manq<-which(is.na(tab))
              
              
            }
            
            tabf<-missForest(tab,maxiter=maxit_method,ntree=ntree)
            tabfinal<-tabf$ximp
            ndm<-length(val_manq)
            diff<-as.matrix(tabfinal)-as.matrix(dtf)
            MSE<-sum(diff^2)/sum(dtf^2)
            list_MSE[i]<-MSE
            list_nbna[i]<-ndm
            
          }
          
        }
        if (method=="EM"){
          for (i in 1:k){
 
            val_manq=rep(0,l)
            
            while(length(val_manq)==0|length(val_manq)==l){
              
              tab<-MCAR(dtf,p)
              val_manq<-which(is.na(tab))
              
            }
            
            tabf<-impute_EM(ds=tab,maxits=maxit_method)
            tabfinal<-tabf
            ndm<-length(val_manq)
            diff<-as.matrix(tabfinal)-as.matrix(dtf)
            MSE<-sum(diff^2)/sum(dtf^2)
            list_MSE[i]<-MSE
            list_nbna[i]<-ndm
            
          }}
        if (method=="SVD"){
          for (i in 1:k){

      
            val_manq=rep(0,l)
            while(length(val_manq)==0|length(val_manq)==l){
             
              tab<-MCAR(dtf,p)
              val_manq<-which(is.na(tab))
            }
            
            tabf<-missSVD(t=tab,n_iter = maxit_method,rg_approx = rg_SVD)
            tabfinal<-tabf
            ndm<-length(val_manq)
            diff<-as.matrix(tabfinal)-as.matrix(dtf)
            MSE<-sum(diff^2)/sum(dtf^2)
            list_MSE[i]<-MSE
            list_nbna[i]<-ndm
            
          }
          
          
        }
        
        
      }
      
      
      
      
      
      
      
      if(miss_data=="MAR"){
      if (method=="miss_forest"){
        for (i in 1:k){
          vect=rep(1,(l+1))
          vmar<-floor(runif(1, min = 1, max = l+1))
          val_manq=rep(0,l)
          
          while(length(val_manq)==0|length(val_manq)==l){
              vect<-rnorm(l+1,0.1,8)
            
            tab<-tableau_MAR(dtf,var_exp =(1:l) ,var_MAR = vmar,dep=vect)
            val_manq<-which(is.na(tab[,vmar]))
            
            
          }
          
          tabf<-missForest(tab,maxiter=maxit_method,ntree=ntree)
          tabfinal<-tabf$ximp
          ndm<-length(val_manq)
          diff<-as.matrix(tabfinal)-as.matrix(dtf)
          MSE<-sum(diff^2)/sum(dtf^2)
          list_MSE[i]<-MSE
          list_nbna[i]<-ndm
          
        }
        
      }
      if(method=="lm"){
        for (i in 1:k){
          vect=rep(1,(l+1))
          vmar<-floor(runif(1, min = 1, max = l+1))
          val_manq=rep(0,l)
          while(length(val_manq)==0|length(val_manq)==l){
            vect<-rnorm(l+1,0.1,8)
            tab<-tableau_MAR(dtf,c(1:l),var_MAR = vmar,dep=vect)
            val_manq<-which(is.na(tab[,vmar]))}
          
          
          ndm<-length(val_manq)
          tabf<-as.data.frame(dtf[-val_manq,])
          colnames(tabf)[vmar]<-"variable_MAR"
          modlm<-lm(variable_MAR~.,data = tabf)
          dt_pred<-dtf[val_manq,]
          dt_pred<-dt_pred[,-vmar]
          
          dt_pred<-cbind(rep(1,ndm),dt_pred)
          
          dt_ver<-dtf[val_manq,]
          dtv<-dt_ver[,vmar]
          coef<-modlm$coefficients
          pred_na<-as.matrix(dt_pred)%*%coef
          res_na<-pred_na-as.matrix(dtv)
          res_na<-res_na^2
          MSE<-sum(res_na)/sum(dtf^2)
          list_MSE[i]<-MSE
          list_nbna[i]<-ndm
          
        }
        
        
      }
      if (method=="EM"){
        for (i in 1:k){
          vect=rep(1,(l+1))
          vmar<-floor(runif(1, min = 1, max = l+1))
          val_manq=rep(0,l)
          
          while(length(val_manq)==0|length(val_manq)==l){
            vect<-rnorm(l+1,0.1,8)
            tab<-tableau_MAR(dtf,c(1:l),var_MAR = vmar,dep=vect)
            val_manq<-which(is.na(tab[,vmar]))
            
            
          }
          
          tabf<-impute_EM(ds=tab,maxits=maxit_method)
          tabfinal<-tabf
          ndm<-length(val_manq)
          diff<-as.matrix(tabfinal)-as.matrix(dtf)
          MSE<-sum(diff^2)/sum(dtf^2)
          list_MSE[i]<-MSE
          list_nbna[i]<-ndm
          
        }}
      if (method=="SVD"){
        for (i in 1:k){
          vect=rep(1,(l+1))
          vmar<-floor(runif(1, min = 1, max = l+1))
          val_manq=rep(0,l)
          while(length(val_manq)==0|length(val_manq)==l){
            vect<-rnorm(l+1,0.1,8)
            tab<-tableau_MAR(dtf,c(1:l),var_MAR = vmar,dep=vect)
            val_manq<-which(is.na(tab[,vmar]))
          }
            
        tabf<-missSVD(t=tab,n_iter = maxit_method,rg_approx = rg_SVD)
        tabfinal<-tabf
        ndm<-length(val_manq)
        diff<-as.matrix(tabfinal)-as.matrix(dtf)
        MSE<-sum(diff^2)/sum(dtf^2)
        list_MSE[i]<-MSE
        list_nbna[i]<-ndm}
        }
      }
      sl_MSE<-mean(list_MSE)
      res<-new(Class="mcdudu",list_mse=list_MSE,list_nna=list_nbna,msetot=sl_MSE)
      return(res)}
  
  
  
  
    porto_MCAR<-MCAR(porto,0.3)
  
    
    
    
    system.time(porto_forest<-missForest(porto_MCAR,ntree=50,maxiter = 10)$ximp)
    porto_alcool_MCAR_reg<-impute_reg(porto_alcool_MCAR,col_na = 12,col_reg = 5)
    porto_alcool_MCAR_reg
    
    
    
    
    
    
    
    ## Calcul de RMSE par monte-carlo répété 100 fois :
    ## SVD itérées
    
    ##Rang 7
    
    Impute_perf_MC(method = "SVD",data.frame(scale(porto[,-1])),rg_SVD =7,maxit_method = 200,
                   k=100,miss_data = "MCAR",p=.3 )@msetot
    ##Rang 9
    
    Impute_perf_MC(method = "SVD",data.frame(scale(porto[,-1])),rg_SVD =9,maxit_method = 200,
                   k=1,miss_data = "MCAR" ,p=.3)@msetot
    
    ##Rang 11
   system.time( MSE_SVD<-Impute_perf_MC(method = "SVD",data.frame(scale(porto[,-1])),rg_SVD =11,maxit_method = 200,
                   k=100,miss_data = "MCAR",p=.3 )@msetot)
    
    
    ## Méthode EM:
    
   system.time(MSE_EM<- Impute_perf_MC(method = "EM",data.frame(scale(porto[,-1])),ntree=50,maxit_method = 100,
                   k=100,miss_data = "MCAR",p=.3 )@msetot)
    
    
    system.time(MSE_reg<-Impute_perf_MC(method = "lm",data.frame(scale(porto[,-1])),
                                        k=100,miss_data = "MCAR" )@msetot)
    
    
    
    
    porto_MCAR<-MCAR(porto,0.3)
    
    system.time(porto_forest<-missForest(porto_MCAR,ntree = 100,maxiter =10)$ximp)
    
    
    Impute_perf_MC(dtf=porto[,c(-1,-14)],method = "miss_forest",k=10,maxit_method = 10,rg_SVD = 2,miss_data = "MCAR",p=0.3)
    porto_MCAR<-MCAR(porto,.3)
    porto_forest<-missForest(porto_MCAR)$ximp
    
    porto_mauvais<-porto
  porto_mauvais<-porto[,porto$quality<=3]
   porto_mauvais<-matrix(0,nrow = nrow(porto),ncol = ncol(porto))
  porto_moyen<-porto_mauvais
  porto_bon<-porto_mauvais
  
  porto$quality1[porto$quality<=3]<-"Mauvais"
  porto$quality1[porto$quality>=7]<-"Bon"
  porto$quality1[porto$quality!="Bon"&porto$quality!="Mauvais"]<-"Moyen"
  porto
 AFD_porto<- desDA(variables = porto[,2:12],group=porto$quality,covar="total")
 
 porto_trou<-MCAR(porto,0.4)
 porto_trou
 missForest(porto_trou)$ximp

 
 
 porto_MCAR[,3]<-MCAR(as.matrix(porto_MCAR[,2]),0.4)
 porto_MCAR[,2]<-porto[,2]
 porto_MCAR
AFD_porto$values
porto_MCAR
    porto_MCAR<-MCAR(porto,0.3)
    porto_MCAR<-cbind(porto[,1:2],porto_MCAR[,4:12],porto$quality)
    porto_reg<-impute_reg(porto_MCAR,col_na = 3:11,col_reg = 2)
    porto_reg
    porto_regression_imputed_AFD<-desDA(porto_reg[,3:11],group=porto$quality,covar = "total")
    porto_regression_imputed_AFD$values
    AFD_porto$values
    
    
  
    help("desDA")
    
    porto
    
  data<-genus
  data
  test<-lm(data$surface~data$altitude)
  
  summary(test)
  
  data[1,]<-colnames(data)
  data<-genus
  
  
  
  
  porto_MCAR<-MCAR(porto[,c(-1,-13)],0.6)
  porto_MCAR
    porto_SVD<-missSVD(porto_MCAR,rg_approx = 10,n_iter = 800)
    porto_SVD
    porto_EM<-impute_EM(porto_MCAR)
    porto_EM
    porto_forest<-missForest(porto_MCAR)
    sum((scale(porto_forest$ximp)-scale(porto[,]))^2)/length(which(is.na(porto_MCAR)))
    
    length(which(is.na(porto_MCAR)))
    
    length(porto_SVD[porto_SVD<=0])
    
 
    
    
   iris_MCAR<-MCAR(iris,0.5) 
    test_SVD<-missSVD(iris_MCAR[,-5])
    
    porto_EM<-impute_reg(porto_MCAR)
  


    
    t<-iris_MCAR[,-5]
    t<-imput_moy(t)
    t-missSVD(iris_MCAR[,-5])
    
    
    
    
    NA_places<-is.na(t)
    t<-imput_moy(t)
    t<-as.matrix(t)
    SVD_obj<- AF(t,nbvrp =min(length(t[,1]),length(t[1,])))
    diag(SVD_obj@lbd)
    t(SVD_obj@U_vects%*%sqrt(diag(SVD_obj@lbd))%*%t(SVD_obj@W_vects))
    
    
    
    
    
    
porto_MCAR<-MCAR(porto[,-c(1,13)])   
 porto_EM<-impute_EM( porto_MCAR)
  
 porto_comp_EM<-cbind(porto[,1],porto_EM,porto[,13])  
 colnames(porto_comp_EM)<-  colnames(porto)
  model_EM<-lm(data = porto_comp_EM,formula = quality~alcohol  )
   summary( model_EM)
##R2 de 27%
## Les valeurs des tests de students sont    
   
   porto$M
   model_porto<-lm(data=porto$alcool_MCAR,formula=quality~-.)
  summary(model_porto)
  ##R2 de 30% dans le mod?le de base 
  
  porto_MCAR_alcool<-MCAR(as.matrix(porto$alcohol),.3)
  porto_MCAR_alcool<-cbind(porto[,c(-12,-13)],porto_MCAR_alcool,porto$quality)
  colnames(porto_MCAR_alcool)=colnames(porto)
  alcohol_EM<-impute_EM(as.matrix(porto_MCAR_alcool))
  colnames(alcohol_EM)<-"alcohol"
  model_alcool_EM<- lm(formula = porto$quality~alcohol_EM)
  summary(model_alcool_EM)
  model_alcool<-lm(formula = porto$quality~porto$alcohol)
  ##0.197 R2 de base
  porto[,c(-12,-13)]
  dim(porto)
  porto_alcool_imputreg<-imput_reg(porto_MCAR_alcool,col_na = 12,col_reg = 2:11)
  colnames(porto_alcool_imputreg)<-colnames(porto)
  model_alcool_imputreg<-lm(data=porto_alcool_imputreg ,formula = quality~alcohol)
  summary(model_alcool_imputreg)
  ##Le R 2 est pas mal si on prend toutes les variables 0.1863
  
  porto_alcool_MCAR<-MCAR(as.matrix(porto$alcohol))

  porto_alcool_MCAR<-impute_EM(porto_alcool_MCAR)
  porto_alcool_MCAR<-cbind(porto[,c(-12,-13)],porto_alcool_MCAR,porto$quality)
  colnames(porto_alcool_MCAR)<-colnames(porto)
  porto_alcool_MCAR_EM
  model_EM_alcohol<-lm(data=porto_alcool_MCAR_EM, formula = quality~alcohol)
  summary(model_EM_alcohol)
  
  
  ##!! EM est un excelent coup!! le R^2 est meilleur que celui de base, 0.1985
  
  porto_alcool_missForest<-missForest(porto_MCAR_alcool[-1,])
  model_forest<-lm(data=porto_alcool_missForest$ximp,formula = quality~alcohol)
  summary(model_forest)
  ##
Analyse_Factorielle_discriminate<- desDA(porto[,c(-1,-13)],porto$quality,covar = "total")
Analyse_Factorielle_discriminate$power
    Analyse_Factorielle_discriminate_imput?e_par_EM<-desDA(porto_EM[,-12],porto$quality,covar = "total")
  Analyse_Factorielle_discriminate_imput?e_par_EM$power
  Analyse_Factorielle_discriminate_imput?e_par_Regression<-desDA(porto_reg)
      
Test<-MCAR(iris,p=0.2 )
length(which(is.na(Test)))/(150*5)     
    
t<-Test[,-5]  

Test_svd<-as.matrix(Test[,-5])
t
Test_svd
missSVD(t)
is.na(Test[,-5])
imput_moy(Test[,-5])
t
class(Test[,-5])

AF_obj_test<-AF(as.matrix(t))


min(length(t[,1]),length(t[1,]))   
AF(t)      
  class(t)
      
      t_moy<-imput_moy(t)
t_test<-as.matrix(t_moy)
  AF_test<-AF(t_test,nbvrp = min(length(t_test[,1]),length(t_test[1,])))  
AF_test@W_vects%*%diag(AF_test@lbd)%*%t(AF_test@U_vects)




cbind(diag(AF_test@lbd),matrix(0,nrow = length(t[,1]),ncol=abs(length(t[,1])-length(t[1,]))))



sum()


diag(AF_test@lbd)
matrix(0,nrow = length(t[,1]),ncol=abs(length(t[,1])-length(t[1,])))





abs(-1)


dim(t)
    t

    
  

