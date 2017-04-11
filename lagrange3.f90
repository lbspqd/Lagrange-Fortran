program lagrange
    implicit none
    !INTERFACE PARA FUNÇÕES E SUBROTINAS
    interface
        real function f(x)
            implicit none
            real x
        end function

        real function val_pol_lagrange(x,y,xi,n)
            implicit none
            real x(:), y(:), xi, L, Pn
            integer n,i,j
        end function
    end interface
    !VARIÁVEIS UTILIZADAS
    integer, parameter :: n = 12 !Pontos que serão interpolados
    integer, parameter :: ni = 3 !Para os arcos
    integer, parameter :: nn = 61 !Para o gráfico
    real x0, xn, passo, passo2
    real, allocatable :: x(:), y(:), xp(:), yp(:)
    integer i,j, aux, aux2

    ! ALOCANDO MEMÓRIA DOS VETORES
    allocate(x(ni)); allocate(y(ni)); allocate(xp(nn-1)); allocate(yp(nn-1));

    !COMEÇA O PROGRAMA
    x0 = -1.0 !PONTO INICIAL
    xn = 1.0 !PONTO FINAL
    passo = (xn - x0)/float(n-1)

    !CRIANDO O PRIMEIRO VETOR DE INTERPOLAÇÃO
    do i=1,ni
        x(i) = x0 + passo*float(i-1)
        y(i) = f(x(i))
    end do
    passo2 = (xn - x0)/float(nn-1)
    aux = 0
    aux2 = 0
    do i=1,n/ni
        do j=1,(nn-1)/(ni+1)
            xp(j+aux2) = x0 + passo2*float(j-1+aux2)
            yp(j+aux2) = val_pol_lagrange(x,y,xp(j+aux2),ni)
        end do
        aux = aux + 3
        aux2 = aux2 + 15
        do j=1,ni
            x(j) = x0 + passo*float(j-1+aux)
            y(j) = f(x(j))
        end do
    end do

    open(1,file="dados3.txt")
    do i=1,nn-1
        write (1,*) xp(i), yp(i)
    end do
    close(1)
    deallocate(x); deallocate(y); deallocate(xp); deallocate(yp)
end program

real function f(x)
    implicit none
    real x
    f = 1/(1+25*x**2)
end function

real function val_pol_lagrange(x,y,xi,n)
    implicit none
    real x(:), y(:), xi, L, Pn
    integer n,i,j

    Pn = 0
    do i=1,n
        L = 1
        do j=1,n
            if (i .NE. j) then
                L =((xi-x(j))/(x(i) - x(j)))*L ! ALGORITMO PARA CALCULO DOS Li'S
            end if
        end do
        Pn = Pn + y(i)*L
    end do
    val_pol_lagrange = Pn
end function
