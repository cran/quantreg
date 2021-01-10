# Toy fnb routine for multiple taus
subroutine qfnb(n,p,m,a,y,t,r,d,u,wn,wp,B,nit, info)

# Input:
#	n = sample size
#	p = parametric dimension of model
#	m = dimension of tau vector
#	a = p by n design matrix (transposed)
#	y = n dimensional response vector
#	t = m dimensional tau vector
#	r = p dimensional rhs vector
#	d = n dimensional vector of ones
#	u = n dimensional vector of ones
#	wn = n by 9 work array
#	wp = p by p+3 work array
#
# Output:	
#	B = p by m matrix of coefficients

integer n,p,m,nit(3),info
double precision a(p,n), y(n), t(m), B(p,m), r(p)
double precision d(n), u(n), wn(n,9), wp(p,p+3)
double precision zero, one, eps, beta

parameter( zero  =  0.0d0)
parameter( one   =  1.0d0)
parameter( beta  =  0.99995d0)
parameter( eps   =  1.0d-6)
    
do i = 1,m{
	call dgemv('N',p,n,one-t(i),a,p,d,1,zero,r,1)
	call dscal(n,zero,wn,1)
	call daxpy(n,one-t(i),u,1,wn,1)
	call rqfnb(n,p,a,y,r,d,u,beta,eps,wn,wp,nit,info)
	if(info != 0) break
	do j = 1,n{
	    u(j) = one
	    d(j) = one
	}
	call dcopy(p,wp,1,B(1,i),1)
    }
return
end
