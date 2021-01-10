# bootstrap estimation with fnb algorithm and globbing

# Input:
    # n = full sample size
    # p = parametric dimension of model
    # m = number of bootstrap replications
    # a = p by n (transposed) design matrix
    # y = n-vector of responses
    # w = n-vector of weights for bootstrap
    # r = residual n-vector from fit 
    # b = coefficient matrix 
    # band = confidence band (see rq.fit.pfnb)
    # n0 = initial sample size 
    # d = n-vector of ones
    # u = n-vector of ones
    # wn = 9*n work array
    # wp = p(p+3) work array
    # aa = n by p array for globbed design matrix
    # yy = n vector for globbed response matrix
    # slo = n-vector of boolean integers
    # shi = n-vector of boolean integers
    # rhs = p-vector for right hand side 
    # glob = p-vector for glob
    # ghib = p-vector for ghimb
    # nit = 3-vector of iteration counts x R
    # info = flag for convergence x R

subroutine pwxy(n,p,m,a,y,tau,qk,r,b,w,band,n0,d,u,wn,wp,
		aa,yy,slo,shi,rhs,glob,ghib,nit,info)
integer n,p,m,kk(2),nn,n0,nit(5,m),info(m),sumbad,i,j,k,ir
integer loq,hiq, slo(n),shi(n),ifix,ibad
logical notopt
double precision a(p,n),y(n),tau,qk(2),r(n),d(n),u(n),b(p,m),w(n)
double precision wn(n,9), wp(p,(p+3)),band(n)
double precision glob(p),ghib(p),aa(p,n),yy(n),rhs(p)
double precision zero,one,beta,eps,big

parameter(zero = 0.0d0)
parameter(one = 1.0d0)
parameter(beta = 0.99995d0)
parameter(big = 1.0d+10)
parameter(eps = 1.0d-06)

    
do ir = 1,m {
  notopt = .true.
  call grexp(n,w,one)
  nn = n0
  ifix = 0
  ibad = 0
  while (notopt) {
    ibad = ibad + 1
    loq = max0(1, int(n*tau - nn/2.))
    hiq = min0(int(n*tau + nn/2.), n)
    qk(1) = r(loq)
    qk(2) = r(hiq)
    call iphil(n,0,slo)
    call iphil(n,0,shi)
    do i = 1,n{
	if(r(i) < qk(1))
	    slo(i) = 1
	else if(r(i) > qk(2))
	    shi(i) = 1
    }
    while (notopt) {
	ifix = ifix + 1
	call dphil(p,zero,glob)
	call dphil(p,zero,ghib)
	call dphil(n,one,d)
	call dphil(n,one,u)
	k = 0
	do i = 1,n{
	    if(slo(i) == 0 & shi(i) == 0){
		k = k + 1
		call dphil(p,zero,aa(1,k))
		call daxpy(p,w(i),a(1,i),1,aa(1,k),1)
                yy(k) = -y(i)*w(i)
		}
	    else if (slo(i) == 1) {
		do j = 1,p 
		    glob(j) = glob(j) + a(j,i) * w(i)
                }
	    else if (shi(i) == 1) {
		do j = 1,p 
		    ghib(j) = ghib(j) + a(j,i) * w(i)
                }
	}
	call dcopy(p,glob,1,aa(1,k+1),1)
	call dcopy(p,ghib,1,aa(1,k+2),1)
        yy(k+1) =  big
        yy(k+2) = -big
	call dgemv('N',p,k+2,one-tau,aa,p,d,1,zero,rhs,1)
	call dscal(k+2,zero,wn,1)
	call daxpy(k+2,one-tau,u,1,wn,1)
        call rqfnb(k+2,p,aa,yy,rhs,d,u,beta,eps,wn,wp,nit(1,ir),info(ir))
	call dcopy(p,wp,1,b(1,ir),1)
	call dcopy(n,y,1,u,1)
	call dgemv('T',p,n,one,a,p,b(1,ir),1,one,u,1)
	# Check predicted signs of residuals (in u)
	sumbad = 0
	do i = 1,n {
	    if((u(i) > 0) & slo(i) == 1){
		slo(i) = 0
		sumbad = sumbad + 1
	    }
	    if((u(i) < 0) & shi(i) == 1){
		shi(i) = 0
		sumbad = sumbad + 1
	    }
	}
        if(sumbad > 0) {
            if (sumbad > 0.1 * nn) {
                nn = min(2 * nn, n)
                break
            }
        }
        else notopt = .false.
    }
    nit(4,ir) = ifix
    nit(5,ir) = ibad
  }
}
return
end

