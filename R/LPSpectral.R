LPSpectral <-
function(W,k,m=8){
        if(k>m){
            k=m;
            warning("k is greater than m, returning first m approximations.")
        }
        P=W/sum(W)
px <- as.vector(apply(P,2,"sum"))
S <- as.matrix(LP.basis(px,m+1))[,1:m]
L <- t(S)%*%P%*%S
        L.svd <-svd(L)
        T <- S%*%L.svd$u[,1:k]
R<-list()
R$LP   <- L
R$Phi  <- T
        R$sval <- L.svd$d[1:k]
return(R)
}
