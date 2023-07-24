
library("missForest")
library(dplyr)

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
            MSE<-sum(diff^2)/ndm
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
            MSE<-sum(diff^2)/ndm
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
            MSE<-sum(diff^2)/ndm
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
          MSE<-sum(diff^2)/ndm
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
          MSE<-sum(res_na)/ndm
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
          MSE<-sum(diff^2)/ndm
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
        MSE<-sum(diff^2)/ndm
        list_MSE[i]<-MSE
        list_nbna[i]<-ndm}
        }
      }
      sl_MSE<-mean(list_MSE)
      res<-new(Class="mcdudu",list_mse=list_MSE,list_nna=list_nbna,msetot=sl_MSE)
      return(res)}
  
  

    
  

