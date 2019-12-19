#library(datasets)
rm(list=ls())
# creating inity vector
uniVc = c(1.0,1.0,1.0,1.0)
# creating the asset matrix
assetMtx <- data.frame( asset = character(0) ,mu= numeric(0), std= numeric(0))
# populating table
assetMtx <- rbind(assetMtx, data.frame(asset = "A", mu =0.04 ,std =0.07 ))
assetMtx <- rbind(assetMtx, data.frame(asset = "B", mu =0.08 ,std =0.12 ))
assetMtx <- rbind(assetMtx, data.frame(asset = "C", mu =0.12 ,std =0.18 ))
assetMtx <- rbind(assetMtx, data.frame(asset = "D", mu =0.15 ,std =0.26 ))

# creating correlation matrix
corrMtx = matrix(nrow = 4, ncol=4, byrow=TRUE)
# populating correlation matrix
corrMtx[1,]=  c(1.0, 0.2, 0.5,0.3)
corrMtx[2,]=  c(0.2, 1.0, 0.7, 0.4)
corrMtx[3,]=  c(0.5, 0.7, 1.0, 0.9)
corrMtx[4,]=  c(0.3, 0.4, 0.9, 1.0)

# creating standart desviation Matrix
stdMtx = matrix(nrow = 4, ncol=4, byrow=TRUE)
# populating standart desviation Matrix
for(i in 1:nrow(corrMtx))
{
  for(j in 1:ncol(corrMtx))
  {
    if(i==j)
    {
      stdMtx[i,j] = assetMtx[i,3]
    } 
    else 
    {
      stdMtx[i,j] = 0.0
    }
    
  }
}

# creating covariances Matrix
# the covariance matrix is obtained by the formula bellow.
# cor(x,y) = cov(x,y)/(std(x)*std(y))
# cov(x,y) = cor(x,y)*(std(x)*std(y))
covMtx =stdMtx %*%(corrMtx%*%stdMtx)
print("Covariance matrix")
print(covMtx)

##################################################################
##################################################################
####    QUESTION  2.
##################################################################
##################################################################

#   2.a) Obtain analytical solution for the the Lagrangian multiplier and optimal allocations w*
#        See resolution for the analitical solution at notebook sheet A_Optimal_Portfolio_Allocations_Exerc_2.jpg


#   2.b) Consider the following optimization task for a targeted 
#        return m = 10%, for which the net of
#        allocations invested (borrowed) in a risk-free asset:

#calculate w*
#step 1 parametrizing solution
targR=0.03 # target return
riskF=0.005
premVc = matrix(nrow = 4, ncol=1, byrow=TRUE)
premVc = assetMtx[,2] -riskF

#Calculating optimal allocations w*
#step 2 calculating scalar Var_A  
Var_A = t(premVc)%*%solve(covMtx)%*%premVc
#step 3 calculating w*
w_optmal = ((targR -riskF )/Var_A[1,1])*solve(covMtx)%*%(premVc)

#step 2 calculating portfolio standard desviation
varPortf = t(w_optmal)%*%covMtx%*%w_optmal
stdPortf = sqrt(varPortf)
print('w* : ')
print(w_optmal)
print('std dev portfolio : ')
print( stdPortf )
