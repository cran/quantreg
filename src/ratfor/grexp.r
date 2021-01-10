
subroutine grexp(n, x, a) 
integer i,n
double precision x(n),a
call fseedi() 
do i = 1,n
   call frexp(x(i), a) 
call fseedo() 
return
end 


