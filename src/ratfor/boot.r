#rq for multple y's
subroutine rqs(m,n,k,m5,n2,a,b,t,toler,ift,x,e,s,
		wa,wb,nsol,sol,dsol,lsol)
double precision b(m,k),sol(n2,nsol),a(m,n),x(n,k)
double precision wa(m5,n2),wb(m),e(m),dsol(m,nsol)
double precision t,toler
integer m,n,k,m5,n2,ift,nsol,lsol
integer s(m)
do i=1,k 
	call rq1(m,n,m5,n2,a,b(1,i),t,toler,ift,x(1,i),e,s,
		wa,wb,nsol,sol,dsol,lsol)
return
end
#parzen, wei and ying's bootstrap
subroutine pwy(m,n,k,m5,n2,a,c,b,t,toler,ift,x,e,s,
		wa,wb,nsol,sol,dsol,lsol)
double precision b(m),sol(n2,nsol),a(k,n),x(n,k)
double precision wa(m5,n2),wb(m),e(m),dsol(m,nsol),c(m,n)
double precision t,toler
integer m,n,k,m5,n2,ift,nsol,lsol
integer s(m)
do i=1,k{
	call dcopy(n,a(i,1),k,c(m,1),m)
	call rq1(m,n,m5,n2,c,b,t,toler,ift,x(1,i),e,s,wa,wb,nsol,sol,dsol,lsol)
	}
return
end
#ratfor outer loop for xy-pairs rq bootstrap
#notation is horrendous 
#   ratfor   R-function
#______________________
#	m -> n
#	n -> p
#	k -> R
#	mofn -> m
subroutine xys(mofn,m,n,k,mofn5,n2,a,b,t,toler,ift,x,e,s,
	wa,wb,nsol,sol,dsol,lsol,aa,bb,ss)
double precision b(m),sol(n2,nsol),a(m,n),x(n,k),wa(mofn5,n2),wb(mofn),e(mofn),dsol(mofn,nsol)
double precision aa(mofn,n),bb(mofn)
double precision t,toler
integer ss(mofn,k),s(mofn),mofn,mofn5
do i=1,k {
	do ii=1,mofn{
		bb(ii)=b(ss(ii,i))
		do jj=1,n{
			aa(ii,jj)=a(ss(ii,i),jj)
			}
		}
	call rq1(mofn,n,mofn5,n2,aa,bb,t,toler,ift,x(1,i),e,s,wa,wb,nsol,sol,dsol,lsol)
	}
return
end
#does a matrix multiply to make Y matrix for heqf bootstrap
subroutine heqfy(n,p,r,x,b,y)
integer n,p,r
double precision x(n,p),b(p,n,r),y(n,r)
do i=1,r{
	do j=1,n{
		y(j,i)=ddot(p,x(j,1),n,b(1,j,i),1)
		}
	}
return
end
