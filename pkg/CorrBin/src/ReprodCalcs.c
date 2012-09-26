
SEXP ReprodEstimates(SEXP nvec, SEXP rvec, SEXP freqvec)
{
 double *theta, *thetanew, sqerror, **marg;
 int i, maxsize, nr, ntot, r, ri, ni, n, fri, start;
 SEXP res;
 const double eps=1e-8;
 
 nr = LENGTH(nvec);
 maxsize = 0;
 ntot = 0;
 for (i=0; i<nr; i++){
    if (INTEGER(nvec)[i]>maxsize){
       maxsize = INTEGER(nvec)[i];
    }
    ntot += INTEGER(freqvec)[i];
 }
       
 theta = malloc((maxsize+1) * sizeof(double));
 thetanew = malloc((maxsize+1) * sizeof(double));
 //starting values
 for (r=0; r<=maxsize; r++){
    theta[r] = 1.0/(maxsize+1);
 }
 sqerror = 1;
 //EM update
 while (sqerror>eps){
    sqerror = 0;
    marg = Marginals(theta, maxsize);
    for (r=0; r<=maxsize; r++) thetanew[r] = 0;
    for (i=0; i<nr; i++){
       ri = INTEGER(rvec)[i];
       ni = INTEGER(nvec)[i];
       fri = INTEGER(freqvec)[i];
       for (r=ri; r<=maxsize-ni+ri; r++){
          thetanew[r] += choose(ni,ri)*choose(maxsize-ni,r-ri)*theta[r]*fri*1.0/
                         marg[ni][ri];
       }
    }
    for (r=0; r<=maxsize; r++){
       thetanew[r] = thetanew[r]/(ntot*choose(maxsize,r)*1.0);
       sqerror += R_pow_di(thetanew[r]-theta[r],2);
       theta[r] = thetanew[r];
    }
    for(n = 0; n <= maxsize; n++) free(marg[n]);
     free(marg);
 }

 PROTECT(res = allocMatrix(REALSXP, maxsize+1, maxsize));
 marg = Marginals(theta, maxsize);
 for (n=1, start=0; n<=maxsize; n++, start+=(maxsize+1)){
    for (r=0; r<=n; r++){
       REAL(res)[start+r] = marg[n][r];
    }
    for (r=n+1; r<=maxsize; r++){
       REAL(res)[start+r] = NA_REAL;
    }
 }
 
 for(n = 0; n <= maxsize; n++) free(marg[n]);
 free(marg);
 free(theta);
 free(thetanew);
 
 UNPROTECT(1);
 return res;
   
}
