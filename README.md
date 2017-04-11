# Lagrange-Fortran

SUBROUTINE intpol(n,g,xi,yi,x,y,opt)
!
! subrotina para a interpolação polinomial
! ou de Lagrange de grau g usando o algoritmo de Neville
! (Numerical Recipes, Cap. 3)
! INPUT 
! xi,yi -> vectores com n elementos contendo os valores de entrada
! da funcao a interpolar ( yi=f(xi) )
! x -> valor da abcissa do ponto onde quer saber o valor da funcao 
! OUTPUT
! y -> valor de f(x)
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: n,g
DOUBLE PRECISION, DIMENSION(n),INTENT(IN) :: xi,yi
DOUBLE PRECISION,INTENT(IN) :: x
DOUBLE PRECISION,INTENT(OUT) :: y
INTEGER, DIMENSION(1) :: iminloc
DOUBLE PRECISION, DIMENSION(0:g,n) :: c,d
DOUBLE PRECISION, DIMENSION(g+1) :: p

INTEGER :: opt, nx, i, m

! n tem de ser tal que n >= g+1
if (n <= g) then 
  write(*,*) 'o numero de pontos dado eh inferior ao grau do polinomio mais um'
  return
end if

iminloc=minloc(abs(x-xi)) 
nx=iminloc(1)



if ( x-xi(nx) < 0.d0 ) nx=nx-1
if ( nx+g > n) nx=n-g

! se opt /= 0, nx=opt -> permite desenhar os polinómios
if (opt /=0) nx=opt


if ( dabs(xi(nx)-x) <= 1.d-16 ) then
  y=yi(nx)
  return
end if 

!
! inicializacao dos c's e d's - > c(0,i)=d(0,i)=yi(i)
  c(0,:)=yi
  d(0,:)=yi
  p=yi(nx:nx+g)
  
do m=0,g-1
  do i=nx, nx+g-1-m
    c(m+1,i)=(xi(i)-x)*(c(m,i+1)-d(m,i))/(xi(i)-xi(i+m+1))  
    d(m+1,i)=(xi(i+m+1)-x)*(c(m,i+1)-d(m,i))/(xi(i)-xi(i+m+1))  
    p(i-nx+1)=0.5d0*(p(i-nx+1)+p(i-nx+2)+c(m+1,i)+d(m+1,i)) 
  end do
end do 
y=p(1)

END SUBROUTINE intpol
