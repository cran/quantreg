C Output from Public domain Ratfor, version 1.05
      subroutine spwy(n,m,k,b1,wu,nnza,a,ja,ia,ao,jao,iao,nnzdmax,d,jd,i
     *d, dsub,jdsub,nnzemax,e,je,ie,nsubmax,lindx,xlindx, nnzlmax,lnz,xl
     *nz,iw,iwmax,iwork,xsuper,tmpmax, tmpvec, wwm,wwn,cachsz,level,x,s,
     *u,c,y,b,small,ierr,maxit,time)
      integer nnza,m,n,nnzdmax,nnzemax,iwmax, nnzlmax,nsubmax,cachsz,lev
     *el,tmpmax,ierr,maxit,time(7), ja(nnza),jao(nnza),jdsub(nnzemax+1),
     *jd(nnzdmax), ia(n+1),iao(m+1),id(m+1),lindx(nsubmax),xlindx(m+1), 
     *iw(m,5),xlnz(m+1),iwork(iwmax),xsuper(m+1),je(nnzemax),ie(m+1)
      double precision small, a(nnza),ao(nnza),dsub(nnzemax+1),d(nnzdmax
     *), lnz(nnzlmax),c(n),y(m),wwm(m,3),tmpvec(tmpmax), wwn(n,14),x(n),
     *s(n),u(n),e(nnzemax),b(m),b1(m),wu(m,k)
      do23000 i=1,k
      call dcopy(m,wu(1,i),1,a(nnza - m),1)
      call dcopy(m,b1,1,y,1)
      call slpfn(n,m,nnza,a,ja,ia,ao,jao,iao,nnzdmax,d,jd,id, dsub,jdsub
     *,nsubmax,lindx,xlindx,nnzlmax,lnz, xlnz,iw(1,1),iw(1,2),iwmax,iwor
     *k,iw(1,3),iw(1,4), xsuper,iw(1,5),tmpmax,tmpvec,wwm(1,2),cachsz, l
     *evel,x,s,u,c,y,b,wwn(1,1),wwn(1,2),wwn(1,3), wwn(1,4),nnzemax,e,je
     *,ie,wwm(1,3),wwn(1,5),wwn(1,6), wwn(1,7),wwn(1,8),wwn(1,9),wwn(1,1
     *0),wwn(1,11), wwn(1,12),wwn(1,13),wwn(1,14),wwm(1,1),small,ierr,ma
     *xit, time)
      call dcopy(m,y,1,wu(1,i),1)
23000 continue
23001 continue
      return
      end
