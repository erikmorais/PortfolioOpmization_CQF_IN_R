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
####    QUESTION  3.
##################################################################
##################################################################

# THE TANGENCY PORTFOLIO
# the tangency portifolio is the optmal portfolio developed in 2 where there is a new constraint W'1 =1.
# All formulae below were extracted the lecture 2

# step 1: caculate the scalar variables, that are used in the formulae
A = t(uniVc)%*%solve(covMtx)%*%uniVc
B = t(assetMtx$mu)%*%solve(covMtx)%*%uniVc
C = t(assetMtx$mu)%*%solve(covMtx)%*%assetMtx$mu
risk_free = 0.005

# step 2: calculate the vector weight of optmized tangent portfolio where sum of w =100%
W_t_Vec = (solve(covMtx)%*%(assetMtx$mu -risk_free )) / ( B[1,1] -A[1,1]*risk_free) 

# step 3: calculate mean tangent portfolio
# two ways to calculate
meant_t = (C - B*risk_free)/( B -A*risk_free)
print(meant_t)
# or 
meant_t = t(W_t_Vec)%*%assetMtx$mu
print(meant_t)

# step 4: calculate std  portfolio
# two ways to calculate
std_t = sqrt( (C - 2*risk_free*B + (risk_free^2)*A  ) / ( B- risk_free*A )^2 )
print(std_t)
#or 
std_t =  sqrt(t(W_t_Vec)%*%covMtx%*%W_t_Vec)
print(std_t)
