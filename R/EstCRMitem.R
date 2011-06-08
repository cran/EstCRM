EstCRMitem <-
function(data,max.item,min.item,max.EMCycle=500,converge=.01) {  #Start Main function

n=ncol(data)  
N=nrow(data) 

if(is.data.frame(data)==FALSE) stop("The input response data is not a data frame.
Please use as.data.frame() and convert your response data to a data frame object before the analysis")

if(length(max.item)!=length(min.item)) stop("The length of max.item vector is not equal to the length of 
min.item vector. Please check your inputs")

if(dim(data)[2]!=length(max.item)) stop("The number of columns in the data is not equal to the length of max.item vector")

if(dim(data)[2]!=length(min.item)) stop("The number of columns in the data is not equal to the length of min.item vector")

for(i in 1:n) {
if(max(na.omit(data[,i]))> max.item[i]) stop("The column ",i," has values higher than the maximum available score in the
user specified max.item vector. Please check and clean your data.")
}

for(i in 1:n) {
if(min(na.omit(data[,i]))< min.item[i]) stop("The column ",i," has values smaller than the minimum available score in the
user specified min.item vector. Please check and clean your data.")
}

if(is.numeric(data[,i])!=TRUE) stop("The column vectors are not numeric. Please check your data")

desc <- as.data.frame(matrix(nrow=n,ncol=4))
colnames(desc) <- c("Mean","SD","Min","Max")
for(i in 1:n) {
desc[i,1]=mean(data[,i],na.rm=TRUE)
desc[i,2]=sd(data[,i],na.rm=TRUE)
desc[i,3]=min(data[,i],na.rm=TRUE)
desc[i,4]=max(data[,i],na.rm=TRUE)
}

for(i in 1:n){data[,i]= data[,i]-min.item[i]}

max.item <- max.item-min.item
min.item <- c(0,0,0,0)

for(i in 1:n) {
if(length(which(data[,i]==max.item[i]))!=0) {
data[which(data[,i]==max.item[i]),i]=max.item[i]-.01
}
if(length(which(data[,i]==0))!=0) {
data[which(data[,i]==0),i]=.01
}
}
data.original <- data

for(i in 1:n){data[,i]= log(data[,i]/(max.item[i]-data[,i]))}

desc$Mean.of.z <- NA
desc$SD.of.z <- NA
for(i in 1:n) {
desc[i,5]=mean(data[,i],na.rm=TRUE)
desc[i,6]=sd(data[,i],na.rm=TRUE)
}
rownames(desc) <- colnames(data.original) 

loglikelihood <- function(ipar,mu,sigma) { 

first.term <- N*sum(log(ipar[,1])+log(ipar[,3]))
sec.term <- sum(rowSums(t(matrix(ipar[,1]^2,nrow=n,ncol=N))*
((t(matrix(ipar[,2],nrow=n,ncol=N))+(data*t(matrix(ipar[,3],nrow=n,ncol=N)))-matrix(rep(mu,n),nrow=N,ncol=n))^2)+
sigma,na.rm=TRUE),na.rm=TRUE)/2
first.term-sec.term
}

estEM <- function(data,ipar) { #start internal function 2 

sigma = 1/(sum(ipar[,1]^2)+1) 
mu <-sigma*rowSums(t(matrix(ipar[,1]^2,ncol=N,nrow=n))*(t(matrix(ipar[,3],ncol=N,nrow=n))*data+t(matrix(ipar[,2],ncol=N,nrow=n)))) # Equation 20 in Shojima's paper

mumean <- mean(mu,na.rm=TRUE)
muvar <- var(mu,na.rm=TRUE)
zijmeanlist <- as.vector(mean(data,na.rm=TRUE))
zijvarlist <- as.vector(diag(var(data,na.rm=TRUE)))
zijmucovlist <- as.vector(cov(data,mu,use="pairwise.complete.obs"))

gamma <- (muvar+sigma)/zijmucovlist  
beta  <- mumean-gamma*zijmeanlist
alpha <- 1/sqrt(gamma^2*zijvarlist-gamma*zijmucovlist)

ipar <- cbind(alpha,beta,gamma)
colnames(ipar) <- c("a","b","alpha")
rownames(ipar) <- colnames(data)

list(ipar,loglikelihood(ipar,mu,sigma)) 
}

ipar <- t(matrix(nrow=3,ncol=n,c(1,0,1)))  
ipar[,2] <- -mean(data,na.rm=TRUE)
mu <- rep(0,N)
sigma <- 1
oldloglik <- loglikelihood(ipar,mu,sigma)

est1 <- estEM(data,ipar)
ipar <- est1[[1]]   
loglik <- est1[[2]] 

d=loglik-oldloglik  
iter <- 1           

loglikelihoods <- c()
ipars <- vector("list",max.EMCycle)

while(abs(d)>converge && iter < max.EMCycle){
loglikelihoods[iter]=loglik
ipars[[iter]] <- ipar
est1 <- estEM(data,ipar)
ipar <- est1[[1]]
d <- est1[[2]]-loglik
loglik <- est1[[2]]
iter <- iter+1
}

itempar <- vector("list",length(loglikelihoods))
for(i in 1:length(loglikelihoods)) itempar[[i]]=ipars[[i]]

for(i in 1:length(loglikelihoods)) itempar[[i]][,3]=1/ipars[[i]][,3]

maximums <- vector("list",length(loglikelihoods))
start <- t(matrix(nrow=3,ncol=n,c(1,0,1)))  
start[,2] <- -mean(data,na.rm=TRUE)
maximums[[1]]<- cbind(max(abs(itempar[[1]]-start)[,1]),max(abs(itempar[[2]]-start)[,2]),max(abs(itempar[[3]]-start)[,3]))
for(i in 1:(length(loglikelihoods)-1)){
maximums[[i+1]]=cbind(max(abs(itempar[[i+1]]-itempar[[i]])[,1]),
max(abs(itempar[[i+1]]-itempar[[i]])[,2]),
max(abs(itempar[[i+1]]-itempar[[i]])[,3]))
}

name <- c()
for(i in 1:length(loglikelihoods)) name[i]=paste("EMCycle",i,"   Largest Parameter Changes=",
round(maximums[[i]],3)[1]," ",round(maximums[[i]],3)[2]," ",round(maximums[[i]],3)[3],
sep="")
names(itempar) <- name

dif <- abs(loglikelihoods[length(loglikelihoods)]-loglikelihoods[length(loglikelihoods)-1])

out <- list(data=data.original,descriptive=desc,param=itempar[[length(loglikelihoods)]],iterations=itempar,dif=dif)
class(out) <- "CRM"
return(out)
}

