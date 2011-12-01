C Output from Public domain Ratfor, version 1.0
      subroutine rqs(m,n,k,m5,n2,a,b,tau,toler,ift,x,e,s,wa,wb)
      double precision b(m,k),a(m,n),x(n,k),e(m),wa(m5,n2),wb(m)
      double precision tau,toler
      integer s(m),m,n,k,m5,n2,ift(k)
      do23000 i=1,k
      call rq0(m,n,m5,n2,a,b(1,i),tau,toler,ift(i),x(1,i),e,s,wa,wb)
23000 continue
23001 continue
      return
      end
