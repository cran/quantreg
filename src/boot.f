      subroutine rqs(m,n,k,m5,n2,a,b,t,toler,ift,x,e,s,wa,wb,nsol,sol,
&     dsol,lsol)
      double precision b(m,k),sol(n2,nsol),a(m,n),x(n,k)
      double precision wa(m5,n2),wb(m),e(m),dsol(m,nsol)
      double precision t,toler
      integer m,n,k,m5,n2,ift,nsol,lsol
      integer s(m)
      do 23000 i=1,k
      call rq1(m,n,m5,n2,a,b(1,i),t,toler,ift,x(1,i),e,s,wa,wb,nsol,sol,
&     dsol,lsol)
23000 continue
      return
      end
      subroutine pwy(m,n,k,m5,n2,a,c,b,t,toler,ift,x,e,s,wa,wb,nsol,sol,
&     dsol,lsol)
      double precision b(m),sol(n2,nsol),a(k,n),x(n,k)
      double precision wa(m5,n2),wb(m),e(m),dsol(m,nsol),c(m,n)
      double precision t,toler
      integer m,n,k,m5,n2,ift,nsol,lsol
      integer s(m)
      do 23002 i=1,k
      call dcopy(n,a(i,1),k,c(m,1),m)
      call rq1(m,n,m5,n2,c,b,t,toler,ift,x(1,i),e,s,wa,wb,nsol,sol,dsol,
&     lsol)
23002 continue
      return
      end
      subroutine xys(mofn,m,n,k,mofn5,n2,a,b,t,toler,ift,x,e,s,wa,wb,
&     nsol,sol,dsol,lsol,aa,bb,ss)
      double precision b(m),sol(n2,nsol),a(m,n),x(n,k),wa(mofn5,n2),wb(
&     mofn),e(mofn),dsol(mofn,nsol)
      double precision aa(mofn,n),bb(mofn)
      double precision t,toler
      integer ss(mofn,k),s(mofn),mofn,mofn5
      do 23004 i=1,k 
      do 23006 ii=1,mofn
      bb(ii)=b(ss(ii,i))
      do 23008 jj=1,n
      aa(ii,jj)=a(ss(ii,i),jj)
23008 continue
23006 continue
      call rq1(mofn,n,mofn5,n2,aa,bb,t,toler,ift,x(1,i),e,s,wa,wb,nsol,
&     sol,dsol,lsol)
23004 continue
      return
      end
      subroutine heqfy(n,p,r,x,b,y)
      integer n,p,r
      double precision x(n,p),b(p,n,r),y(n,r)
      do 23010 i=1,r
      do 23012 j=1,n
      y(j,i)=ddot(p,x(j,1),n,b(1,j,i),1)
23012 continue
23010 continue
      return
      end
