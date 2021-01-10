# multiple tau estimation with fnb algorithm and globbing

# Input:
    # n = full sample size
    # p = parametric dimension of model
    # m = number of quantiles to be estimated
    # a = p by n (transposed) design matrix
    # y = n-vector of responses
    # q = m-vector of quantiles to be estimated
    # r = residual vector from fit for q[1]
    # b = coefficient matrix 
    # band = confidence band (see rq.fit.pfnb)
    # m0 = initial sample size default is n^(2/3) * p^(1/2)
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
    # ghib = p-vector for ghib
    # nit = 3-vector of iteration counts
    # info = flag for convergence
    # Maybe there should be a control object 
    # Possible control parameter m0
    # and/or work storage could be better rationalized.

subroutine pfnb(n,p,m,a,y,q,r,b,band,m0,d,u,wn,wp,
		aa,yy,slo,shi,rhs,glob,ghib,nit,info)
integer n,p,m,kk(2),mm,m0,nit(5,m),info,sumbad
integer i,j,k,slo(n),shi(n),ifix,ibad
logical notopt
double precision a(p,n),y(n),r(n),q(m),d(n),u(n),b(p,m)
double precision wn(n,9), wp(p,(p+3)),band(n),qk(2)
double precision glob(p),ghib(p),aa(p,n),yy(n),rhs(p)
double precision zero,one,tau,beta,eps,big,fm,fn

parameter(zero = 0.0d0)
parameter(one = 1.0d0)
parameter(beta = 0.99995d0)
parameter(big = 1.0d+10)
parameter(eps = 1.0d-06)

    
do iq = 1,m {
  notopt = .true.
  tau = q(iq)
  mm = m0
  ifix = 0
  ibad = 0
  while (notopt) {
    ibad = ibad + 1
    fm = mm
    fn = dble(n)
    kk(1) = int(n * dmax1(1./fn, tau - fm/(2 * fn))) + 1
    kk(2) = int(n * dmin1(tau + fm/(2. * fn), (fn - 1)/fn))
    do i = 1,n
	u(i) = r(i)/band(i)
    call kuantiles(kk,2,n,u)
    qk(1) = u(kk(1))
    qk(2) = u(kk(2))
    call iphil(n,0,slo)
    call iphil(n,0,shi)
    do i = 1,n{
	if(r(i) < (band(i) * qk(1)))
	    slo(i) = 1
	else if(r(i) > (band(i) * qk(2)))
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
		call dcopy(p,a(1,i),1,aa(1,k),1)
                yy(k) = -y(i)
		}
	    else if (slo(i) == 1) {
		do j = 1,p 
		    glob(j) = glob(j) + a(j,i) 
                }
	    else if (shi(i) == 1) {
		do j = 1,p 
		    ghib(j) = ghib(j) + a(j,i) 
                }
	}
	call dcopy(p,glob,1,aa(1,k+1),1)
	call dcopy(p,ghib,1,aa(1,k+2),1)
        yy(k+1) =  big
        yy(k+2) = -big
	call dgemv('N',p,k+2,one-tau,aa,p,d,1,zero,rhs,1)
	call dscal(k+2,zero,wn,1)
	call daxpy(k+2,one-tau,u,1,wn,1)
        call rqfnb(k+2,p,aa,yy,rhs,d,u,beta,eps,wn,wp,nit(1,iq),info)
	call dcopy(p,wp,1,b(1,iq),1)
	call dcopy(n,y,1,r,1)
	call dgemv('T',p,n,one,a,p,b(1,iq),1,one,r,1)
	# Check predicted signs of residuals (in r)
	sumbad = 0
	do i = 1,n {
	    if((r(i) > 0) & slo(i) == 1){
		slo(i) = 0
		sumbad = sumbad + 1
	    }
	    if((r(i) < 0) & shi(i) == 1){
		shi(i) = 0
		sumbad = sumbad + 1
	    }
	}
        if(sumbad > 0) {
            if (sumbad > 0.1 * mm) {
                mm = min(2 * mm, n)
                break
            }
        }
        else notopt = .false.
    }
    nit(4,iq) = ifix
    nit(5,iq) = ibad
  }
}
return
end
subroutine iphil(n,a,v)
integer n,i,a,v(n)
do i = 1,n
    v(i) = a
return
end
subroutine dphil(n,a,v)
integer i,n
double precision a,v(n)
do i = 1,n
    v(i) = a
return
end
