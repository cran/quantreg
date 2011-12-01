      subroutine rq(m,nn,m5,n2,n4,a,b,t,toler,ift,x,e,s,wa,wb,nsol,
&     ndsol,sol,dsol,lsol,h,qn,cutoff,ci,tnmat,big,lci1)
      integer i,j,k,kl,kount,kr,l,lsol,m,m1,m2,m3,m4,m5,ift
      integer n,n1,n2,nsol,ndsol,out,s(m),h(nn,nsol)
      integer nn,n4,idxcf
      logical stage,test,init,iend,lup
      logical lci1,lci2,skip
      real a1,aux,b1,big,d,dif,pivot,smax,t,t0,t1,tnt
      real min,max,toler,zero,half,one,two
      real b(m),sol(n2,nsol),a(m,nn),x(nn),wa(m5,n4),wb(m)
      real sum,e(m),dsol(m,ndsol)
      real qn(nn),cutoff,ci(4,nn),tnmat(4,nn),tnew,told,tn
      data zero/0.00/
      data half/0.50/
      data one/1.00/
      data two/2.00/
      n=nn
      ift = 0
      wa(m+2,nn+1) = one
      if(.not.(m5.ne.m+5))goto 23000
      ift = 3
23000 continue
      if(.not.(n2.ne.n+2))goto 23002
      ift = 4
23002 continue
      if(.not.(n4.ne.n+4))goto 23004
      ift = 5
23004 continue
      if(.not.(m.le.zero.or.n.le.zero))goto 23006
      ift = 6
23006 continue
      if(.not.(ift.le.two))goto 23008
      iend = .true.
      lci2 = .false.
      lup = .true.
      skip = .false.
      idxcf = 0
      tnew = zero
      tn = zero
      m1 = m+1
      n1 = n+1
      n3 = n+3
      m2 = m+2
      m3 = m+3
      m4 = m+4
      do 23010 j = 1,n
      x(j) = zero
23010 continue
      do 23012 i = 1,m
      e(i) = zero
23012 continue
      if(.not.(t.lt.zero.or.t.gt.one))goto 23014
      t0 = one/float(m)-toler
      t1 = one-t0
      t = t0
      iend = .false.
      lci1 = .false.
23014 continue
23016 continue
      do 23019 i = 1,m 
      k = 1
      do 23021 j = 1,nn
      if(.not.(k.le.nn))goto 23023
      if(.not.(j.eq.idxcf))goto 23025
      skip = .true.
      goto 23026
23025 continue
      skip = .false.
23026 continue
      if(.not.(.not.skip))goto 23027
      wa(i,k) = a(i,j)
      k = k+1
23027 continue
23023 continue
23021 continue
      wa(i,n4) = n+i
      wa(i,n2) = b(i)
      if(.not.(idxcf .ne. 0))goto 23029
      wa(i,n3) = tnew*a(i,idxcf)
      goto 23030
23029 continue
      wa(i,n3) = zero
23030 continue
      wa(i,n1) = wa(i,n2)-wa(i,n3)
      if(.not.(wa(i,n1).lt.zero))goto 23031
      do 23033 j = 1,n4
      wa(i,j) = -wa(i,j)
23033 continue
23031 continue
23019 continue
      do 23035 j = 1,n 
      wa(m4,j) = j
      wa(m2,j) = zero
      wa(m3,j) = zero
      do 23037 i = 1,m 
      aux = sign(one,wa(m4,j))*wa(i,j)
      wa(m2,j) = wa(m2,j)+aux*(one-sign(one,wa(i,n4)))
      wa(m3,j) = wa(m3,j)+aux*sign(one,wa(i,n4))
23037 continue
      wa(m3,j) = two*wa(m3,j)
23035 continue
      dif = zero
      init = .false.
      if(.not.(.not.lci2))goto 23039
      do 23041 k = 1,n 
      wa(m5,k) = zero
      do 23043 i = 1,m
      wa(m5,k) = wa(m5,k)+a(i,k)
23043 continue
      wa(m5,k) = wa(m5,k)/float(m)
23041 continue
23039 continue
      lsol = 1
      kount = 0
23045 continue
      do 23048 j = 1,n
      wa(m1,j) = wa(m2,j)+wa(m3,j)*t
23048 continue
      if(.not.(.not.init))goto 23050
      stage = .true.
      kr = 1
      kl = 1
      go to 30
23050 continue
23052 continue
      stage = .false.
23055 continue
      max = -big
      do 23058 j = kr,n 
      d = wa(m1,j)
      if(.not.(d.lt.zero))goto 23060
      if(.not.(d.gt.(-two)))goto 23062
      goto 23058
23062 continue
      d = -d-two
23060 continue
      if(.not.(d.gt.max))goto 23064
      max = d
      in = j
23064 continue
23058 continue
      if(.not.(max.le.toler))goto 23066
      goto 23054
23066 continue
      if(.not.(wa(m1,in).le.zero))goto 23068
      do 23070 i = 1,m4
      wa(i,in) = -wa(i,in)
23070 continue
      wa(m1,in) = wa(m1,in)-two
      wa(m2,in) = wa(m2,in)-two
23068 continue
23072 continue
      k = 0
      do 23075 i = kl,m 
      d = wa(i,in)
      if(.not.(d.gt.toler))goto 23077
      k = k+1
      wb(k) = wa(i,n1)/d
      s(k) = i
      test = .true.
23077 continue
23075 continue
23079 continue
      if(.not.(k.le.0))goto 23082
      test = .false.
      goto 23083
23082 continue
      min = big
      do 23084 i = 1,k
      if(.not.(wb(i).lt.min))goto 23086
      j = i
      min = wb(i)
      out = s(i)
23086 continue
23084 continue
      wb(j) = wb(k)
      s(j) = s(k)
      k = k-1
23083 continue
      if(.not.(.not.test.and.stage))goto 23088
      goto 23081
23088 continue
      if(.not.(.not.test))goto 23090
      goto 23047
23090 continue
      pivot = wa(out,in)
      if(.not.(wa(m1,in)-pivot-pivot.le.toler))goto 23092
      go to 10
23092 continue
      do 23094 j = kr,n3 
      d = wa(out,j)
      wa(m1,j) = wa(m1,j)-d-d
      wa(m2,j) = wa(m2,j)-d-d
      wa(out,j) = -d
23094 continue
      wa(out,n4) = -wa(out,n4)
23080 goto 23079
23081 continue
      do 23096 i = 1,m4 
      d = wa(i,kr)
      wa(i,kr) = wa(i,in)
      wa(i,in) = d
23096 continue
      kr = kr+1
      go to 20
10    do 23098 j = kr,n3
      if(.not.(j.ne.in))goto 23100
      wa(out,j) = wa(out,j)/pivot
23100 continue
23098 continue
      do 23102 i = 1,m3
      if(.not.(i.ne.out))goto 23104
      d = wa(i,in)
      do 23106 j = kr,n3
      if(.not.(j.ne.in))goto 23108
      wa(i,j) = wa(i,j)-d*wa(out,j)
23108 continue
23106 continue
23104 continue
23102 continue
      do 23110 i = 1,m3
      if(.not.(i.ne.out))goto 23112
      wa(i,in) = -wa(i,in)/pivot
23112 continue
23110 continue
      wa(out,in) = one/pivot
      d = wa(out,n4)
      wa(out,n4) = wa(m4,in)
      wa(m4,in) = d
      kount = kount+1
      if(.not.(.not.stage))goto 23114
      goto 23074
23114 continue
      kl = kl+1
      do 23116 j = kr,n4 
      d = wa(out,j)
      wa(out,j) = wa(kount,j)
      wa(kount,j) = d
23116 continue
20    if(.not.(kount+kr.eq.n1))goto 23118
      goto 23057
23118 continue
30    max = -one
      do 23120 j = kr,n
      if(.not.(abs(wa(m4,j)).le.n))goto 23122
      d = abs(wa(m1,j))
      if(.not.(d.gt.max))goto 23124
      max = d
      in = j
23124 continue
23122 continue
23120 continue
      if(.not.(wa(m1,in).lt.zero))goto 23126
      do 23128 i = 1,m4
      wa(i,in) = -wa(i,in)
23128 continue
23126 continue
23073 goto 23072
23074 continue
23056 goto 23055
23057 continue
23053 goto 23052
23054 continue
      if(.not.(kr.eq.1))goto 23130
      do 23132 j = 1,n 
      d = abs(wa(m1,j))
      if(.not.(d.le.toler.or.two-d.le.toler))goto 23134
      ift = 1
      wa(m2,nn+1) = zero
      go to 80
23134 continue
23132 continue
23130 continue
80    kount = 0
      sum = zero
      if(.not.(.not.lci2))goto 23136
      do 23138 i = 1,kl-1 
      k = wa(i,n4)*sign(one,wa(i,n4))
      x(k) = wa(i,n1)*sign(one,wa(i,n4))
23138 continue
23136 continue
      do 23140 i = 1,n 
      kd = abs(wa(m4,i))-n
      dsol(kd,lsol) = one+wa(m1,i)/two
      if(.not.(wa(m4,i).lt.zero))goto 23142
      dsol(kd,lsol) = one-dsol(kd,lsol)
23142 continue
      if(.not.(.not.lci2))goto 23144
      sum = sum+x(i)*wa(m5,i)
      sol(i+2,lsol) = x(i)
      h(i,lsol) = kd
23144 continue
23140 continue
      do 23146 i = kl,m 
      kd = abs(wa(i,n4))-n
      if(.not.(wa(i,n4).lt.zero))goto 23148
      dsol(kd,lsol) = zero
23148 continue
      if(.not.(wa(i,n4).gt.zero))goto 23150
      dsol(kd,lsol) = one
23150 continue
23146 continue
      if(.not.(.not.lci2))goto 23152
      sol(1,lsol) = smax
      sol(2,lsol) = sum
      do 23154 i=1,m
      dsol(i,lsol+1) = dsol(i,lsol)
23154 continue
23152 continue
      if(.not.(lci2))goto 23156
      a1 = zero
      do 23158 i = 1,m 
      a1 = a1+a(i,idxcf)*(dsol(i,lsol)+t-one)
23158 continue
      tn = a1/sqrt(qn(idxcf)*t*(one-t))
      if(.not.(abs(tn).lt.cutoff))goto 23160
      if(.not.(lup))goto 23162
      smax = big
      goto 23163
23162 continue
      smax = -big
23163 continue
      do 23164 i =1,kl-1 
      k = wa(i,n4)*sign(one,wa(i,n4))
      sol(k,1) = wa(i,n2)*sign(one,wa(i,n4))
      sol(k,2) = wa(i,n3)*sign(one,wa(i,n4))/tnew
23164 continue
      do 23166 i = kl,m 
      a1 = zero
      b1 = zero
      k = wa(i,n4)*sign(one,wa(i,n4))-n
      l = 1
      do 23168 j = 1,n
      if(.not.(j.eq.idxcf))goto 23170
      l = l+1
23170 continue
      a1 = a1 + a(k,l)*sol(j,1)
      b1 = b1 + a(k,l)*sol(j,2)
      l = l+1
23168 continue
      tnt = (b(k)-a1)/(a(k,idxcf)-b1)
      if(.not.(lup))goto 23172
      if(.not.(tnt.gt.tnew))goto 23174
      if(.not.(tnt.lt.smax))goto 23176
      smax = tnt
      out = i
23176 continue
23174 continue
      goto 23173
23172 continue
      if(.not.(tnt.lt.tnew))goto 23178
      if(.not.(tnt.gt.smax))goto 23180
      smax = tnt
      out = i
23180 continue
23178 continue
23173 continue
23166 continue
      if(.not.(lup))goto 23182
      told = tnew
      tnew = smax + toler
      ci(3,idxcf) = told - toler
      tnmat(3,idxcf) = tn
      if(.not.(.not.(tnew .lt. big-toler)))goto 23184
      ci(3,idxcf) = big
      ci(4,idxcf) = big
      tnmat(3,idxcf) = tn
      tnmat(4,idxcf) = tn
      lup = .false.
      go to 70
23184 continue
      goto 23183
23182 continue
      told = tnew
      tnew = smax - toler
      ci(2,idxcf) = told + toler
      tnmat(2,idxcf) = tn
      if(.not.(.not.(tnew .gt. -big+toler)))goto 23186
      ci(2,idxcf) = -big
      ci(1,idxcf) = -big
      tnmat(2,idxcf) = tn
      tnmat(1,idxcf) = tn
      lup = .true.
      go to 60
23186 continue
23183 continue
      do 23188 i = 1,m
      wa(i,n3) = wa(i,n3)/told*tnew
      wa(i,n1) = wa(i,n2) - wa(i,n3)
23188 continue
      do 23190 j = kr,n3
      d = wa(out,j)
      wa(m1,j) = wa(m1,j) -d -d
      wa(m2,j) = wa(m2,j) -d -d
      wa(out,j) = -d
23190 continue
      wa(out,n4) = -wa(out,n4)
      init = .true.
      goto 23161
23160 continue
      if(.not.(lup))goto 23192
      ci(4,idxcf) = tnew - toler
      tnmat(4,idxcf) = tn
      lup = .false.
      go to 70
      goto 23193
23192 continue
      ci(1,idxcf) = tnew + toler
      tnmat(1,idxcf) = tn
      lup = .true.
      go to 60
23193 continue
23161 continue
23156 continue
      if(.not.((iend).and.(.not.lci2)))goto 23194
      go to 40
23194 continue
      if(.not.(.not.lci2))goto 23196
      init = .true.
      lsol = lsol+1
      do 23198 i = 1,m
      s(i) = zero
23198 continue
      do 23200 j = 1,n
      x(j) = zero
23200 continue
      smax = two
      do 23202 j = 1,n 
      b1 = wa(m3,j)
      a1 = (-two-wa(m2,j))/b1
      b1 = -wa(m2,j)/b1
      if(.not.(a1.ge.t))goto 23204
      if(.not.(a1.lt.smax))goto 23206
      smax = a1
      dif = (b1-a1)/two
23206 continue
23204 continue
      if(.not.(b1.gt.t))goto 23208
      if(.not.(b1.lt.smax))goto 23210
      smax = b1
      dif = (b1-a1)/two
23210 continue
23208 continue
23202 continue
      tnt = smax+toler*(one+abs(dif))
      if(.not.(tnt.ge.t1+toler))goto 23212
      iend = .true.
23212 continue
      t = tnt
      if(.not.(iend))goto 23214
      t = t1
23214 continue
23196 continue
23046 goto 23045
23047 continue
      wa(m2,nn+1) = two
      ift = 2
      go to 50
40    if(.not.(lsol.gt.2))goto 23216
      sol(1,1) = zero
      sol(1,lsol) = one
      do 23218 i = 1,m 
      dsol(i,1) = one
      dsol(i,lsol) = zero
      dsol(i,lsol+1) = zero
23218 continue
23216 continue
      l = kl-1
      do 23220 i = 1,l
      if(.not.(wa(i,n1).lt.zero))goto 23222
      do 23224 j = kr,n4
      wa(i,j) = -wa(i,j)
23224 continue
23222 continue
23220 continue
50    sum = zero
      if(.not.(.not.lci2))goto 23226
      do 23228 i = kl,m 
      k = wa(i,n4)*sign(one,wa(i,n4))
      d = wa(i,n1)*sign(one,wa(i,n4))
      sum = sum+d*sign(one,d)*(half+sign(one,d)*(t-half))
      k = k-n
      e(k) = d
23228 continue
      wa(m2,n2) = kount
      wa(m1,n2) = n1-kr
      wa(m1,n1) = sum
23226 continue
      if(.not.(wa(m2,nn+1).eq.two))goto 23230
      goto 23018
23230 continue
      if(.not.(.not.lci1))goto 23232
      goto 23018
23232 continue
      if(.not.(.not.lci2))goto 23234
      lci2 = .true.
      n = nn-1
      n1 = n+1
      n2 = n+2
      n3 = n+3
      n4 = n+4
60    idxcf = idxcf+1
      if(.not.(idxcf.gt.nn))goto 23236
      goto 23018
23236 continue
70    if(.not.(lup))goto 23238
      tnew = x(idxcf)+toler
      told = tnew
      ci(3,idxcf) = x(idxcf)
      tnmat(3,idxcf) = zero
      goto 23239
23238 continue
      tnew = x(idxcf)-toler
      told = tnew
      ci(2,idxcf) = x(idxcf)
      tnmat(2,idxcf) = zero
23239 continue
23234 continue
23017 goto 23016
23018 continue
      do 23240 i=1,m
      dsol(i,lsol) = dsol(i,lsol+1)
23240 continue
23008 continue
      return
      end
