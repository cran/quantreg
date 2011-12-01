      subroutine dsel05( k, n, x)
      integer            k, n
      double precision   x(n)
c
c     Selects the smallest k elements of the array x[1:n].
c     The input array is permuted so that the smallest k elements of
c     x are x(i), i = 1,...,k, (in arbitrary order) and x(k) is the
c     kth smallest element.
c
c     This is a Fortran 77 version of the Algol 68 procedure from
c
c        R.W. Floyd and R.L. Rivest: "Algorithm 489: The Algorithm
c        SELECT---for Finding the $i$th Smallest of $n$ Elements",
c        Comm. ACM 18, 3 (1975) 173,
c
c     including some modifications suggested in
c
c        T. Brown: "Remark on Algorithm 489", ACM Trans. Math.
c        Software 3, 2 (1976), 301-304.
c
c     Array stack(2,nstack) permits up to nstack levels of recursion.
c     For standard parameters cs <= 1 and cutoff >= 600,
c     nstack = 5 suffices for n up to 2**31-1 (maximum integer*4).
      integer            nstack
      parameter         (nstack=10)
      integer            stack(2,nstack)
c
c     Parameters cutoff, cs and csd are as originally proposed.
      integer            cutoff
      parameter         (cutoff=600)
      double precision   cs, csd
      parameter         (cs=0.5d0, csd=0.5d0)
c     Brown's version
c     parameter         (cs=0.5d0, csd=0.1d0)
c
c     Subprograms called:
      intrinsic          dble, exp, log, max, min, sign
c
c     Written by K.C. Kiwiel, 8 March 2006, kiwiel@ibspan.waw.pl
c
c     Local variables:
      integer            i, j, jstack, l, m, r, s, sd
      double precision   dm, swap, v, z
      l=1
      r=n
      jstack=0
c     entry to SELECT( x, n, l, r, k)
c     SELECT will rearrange the values of the array segment x[l:r] so
c     that x(k) (for some given k; l <= k <= r) will contain the
c     (k-l+1)-th smallest value, l <= i <= k will imply x(i) <= x(k),
c     and k <= i <= r will imply x(k) <= x(i).
c     while r > l do
    1 continue
      if (l.ge.r) goto 6
c        The additional test below prevents stack overflow.
         if (r-l.gt.cutoff .and. jstack.lt.nstack) then
c           Use SELECT recursively on a sample of size s to get an
c           estimate for the (k-l+1)-th smallest element into x(k),
c           biased slightly so that the (k-l+1)-th element is
c           expected to lie in the smaller set after partitioning.
            m=r-l+1
            i=k-l+1
            dm=m
            z=log(dm)
            s=cs*exp(2*z/3)+0.5d0
            sd=csd*sqrt(z*s*(1-s/dm))*sign(1d0,i-dm/2)+0.5d0
            if (i.eq.m/2) sd=0
c           Brown's modification
c           sd=csd*sqrt(z*s*(1-s/dm))*(2*i/dm-1)+0.5d0
c           Push the current l and r on the stack.
            jstack=jstack+1
            stack(1,jstack)=l
            stack(2,jstack)=r
c           Find new l and r for the next recursion.
            l=max(dble(l),k-i*(s/dm)+sd)+0.5d0
            r=min(dble(r),k-i*(s/dm)+sd+s)+0.5d0
c           call SELECT( x, n, l, r, k)
            goto 1
         endif
    2    continue
c        Partition x[l:r] about the pivot v := x(k).
         v=x(k)
c        Initialize pointers for partitioning.
         i=l
         j=r
c        Swap x(l) and x(k).
         x(k)=x(l)
         x(l)=v
         if (v.lt.x(r)) then
c           Swap x(l) and x(r).
            x(l)=x(r)
            x(r)=v
         endif
c        while i < j do
    3    continue
         if (i.lt.j) then
c           Swap x(i) and x(j).
            swap=x(j)
            x(j)=x(i)
            x(i)=swap
            i=i+1
            j=j-1
c           Scan up to find element >= v.
    4       continue
            if (x(i).lt.v) then
               i=i+1
               goto 4
            endif
c           Scan down to find element <= v.
    5       continue
            if (x(j).gt.v) then
               j=j-1
               goto 5
            endif
            goto 3
c           end of while i < j do
         endif
         if (x(l).eq.v) then
c           Swap x(l) and x(j).
            swap=x(j)
            x(j)=v
            x(l)=swap
         else
            j=j+1
c           Swap x(j) and x(r).
            swap=x(j)
            x(j)=x(r)
            x(r)=swap
         endif
c        Now adjust l, r so that they surround the subset containing
c        the (k-l+1)-th smallest element.
         if (j.le.k) l=j+1
         if (k.le.j) r=j-1
         goto 1
c        end of while r > l do
    6 continue
c     Exit if the stack is empty.
      if (jstack.eq.0) return
c     Pop l and r from the stack.
      l=stack(1,jstack)
      r=stack(2,jstack)
      jstack=jstack-1
c     Continue as if after a return from a recursive call.
      goto 2
      end
