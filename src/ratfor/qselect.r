#function to compute qth quantile of a sample of n observations
subroutine qselect(n,x,q)
integer n,k,l,r
double precision x(n),q
k=nint(q*n)
l=1
r=n
call select(n,x,l,r,k)
q=x(k)
return
end
#This is a ratfor implementation of the floyd-rivest algorithm--SELECT
#Reference:  CACM 1975, alg #489, p173, algol-68 version
#Translation by Roger Koenker August, 1996.
#As originally proposed mmax=600, and cs=cd=.5
#Calls blas routine dswap
recursive subroutine select(n,x,l,r,k)
integer n,m,l,r,k,ll,rr,i,j,mmax
double precision x(n),z,s,d,t,fm,cs,cd
parameter(cs = 0.5d0)
parameter(cd = 0.5d0)
parameter(mmax = 600)
while(r>l){
	if(r-l>mmax){
		m=r-l+1
		i=k-l+1
		fm = dble(m)
		z=log(fm)
		s=cs*exp(2*z/3)
		d=cd*sqrt(z*s*(m-s)/fm)*sign(1,i-m/2)
		ll=max(l,nint(k-i*s/fm + d))
		rr=min(r,nint(k+(m-i)*s/fm + d))
		call select(n,x,ll,rr,k)
		}
	t=x(k)
	i=l
	j=r
	call dswap(1,x(l),1,x(k),1)
	if(x(r)>t)call dswap(1,x(r),1,x(l),1)
	while(i<j){
		call dswap(1,x(i),1,x(j),1)
		i=i+1
		j=j-1
		while(x(i)<t)i=i+1
		while(x(j)>t)j=j-1
		}
	if(x(l)==t)
		call dswap(1,x(l),1,x(j),1)
	else{
		j=j+1
		call dswap(1,x(j),1,x(r),1)
		}
	if(j<=k)l=j+1
	if(k<=j)r=j-1
	}
return
end
