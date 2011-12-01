C Output from Public domain Ratfor, version 1.0
      subroutine pwy(m,n,k,m5,n2,a,c,b,t,toler,ift,x,e,s, wa,wb)
      double precision b(m),a(k,n),x(n,k)
      double precision wa(m5,n2),wb(m),e(m),c(m,n)
      double precision t,toler
      integer m,n,k,m5,n2,ift
      integer s(m)
      do23000 i=1,k
      call dcopy(n,a(i,1),k,c(m,1),m)
      call rq0(m,n,m5,n2,c,b,t,toler,ift,x(1,i),e,s,wa,wb)
23000 continue
23001 continue
      return
      end
      subroutine xys(mofn,m,n,k,mofn5,n2,a,b,tau,toler,ift,x,e,s, wa,wb,
     *aa,bb,ss)
      double precision b(m),a(m,n),x(n,k)
      double precision wa(mofn5,n2),wb(mofn)
      double precision aa(mofn,n),bb(mofn),e(mofn)
      double precision tau,toler
      integer ss(mofn,k),s(mofn),mofn,m,n,k,mofn5,n2,ift(k)
      do23002 i=1,k 
      do23004 ii=1,mofn
      bb(ii)=b(ss(ii,i))
      do23006 jj=1,n
      aa(ii,jj)=a(ss(ii,i),jj)
23006 continue
23007 continue
23004 continue
23005 continue
      call rq0(mofn,n,mofn5,n2,aa,bb,tau,toler,ift(i),x(1,i),e,s,wa,wb)
23002 continue
23003 continue
      return
      end
      subroutine wxy(m,n,k,m5,n2,a,b,tau,toler,ift,x,e,s,wa,wb,aa,bb,w)
      double precision b(m),a(m,n),x(n,k)
      double precision w(m,k),wa(m5,n2),wb(m)
      double precision aa(m,n),bb(m),e(m)
      double precision tau,toler
      integer s(m),m,n,k,m5,n2,ift(k)
      do23008 i=1,k 
      do23010 ii=1,m
      bb(ii)=b(ii)*w(ii,i)
      do23012 jj=1,n
      aa(ii,jj)=a(ii,jj)*w(ii,i)
23012 continue
23013 continue
23010 continue
23011 continue
      call rq0(m,n,m5,n2,aa,bb,tau,toler,ift(i),x(1,i),e,s,wa,wb)
23008 continue
23009 continue
      return
      end
      subroutine heqfy(n,p,r,x,b,y)
      integer n,p,r
      double precision x(n,p),b(p,n,r),y(n,r)
      do23014 i=1,r
      do23016 j=1,n
      y(j,i)=ddot(p,x(j,1),n,b(1,j,i),1)
23016 continue
23017 continue
23014 continue
23015 continue
      return
      end
