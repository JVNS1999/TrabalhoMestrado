!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                     Programa Cálculo do Monte Carlo                    !         
!                     Programador: João Vitor Nunes                      !                 
!                        Orientador: Lucas Mól                           !               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Module Variaveis
    Implicit None
    Integer :: Ns, Nterma
    Integer :: Nmc
    Integer, Dimension(:), allocatable :: S
    Integer(4) :: cont
    Real(8), Dimension(:,:), allocatable :: Aij
    Real(8), Dimension(:), allocatable :: Bi, Ei,ex,ey
    Real(8) :: beta , Etot, Mx, My, M
End Module Variaveis

Subroutine LerAij()
    Use Variaveis
    Implicit None
    Integer :: I , J
    real(8)::x,y
    Open(20,file = 'Aij.dat', status = 'Old')
    Read(20,*)Ns
    Allocate(Aij(Ns,Ns),ex(Ns),ey(Ns))
    Do i = 1,Ns
        Do j = 1,Ns
            Read(20,*) Aij(i,j)
        End Do
    End Do
    Close(20)  


    Open(10,file = 'resultados.csv')
    Read(10,*)
    Do i = 1,Ns
        Read(10,*)x,y,ex(i),ey(i)
    End Do
    Close(10)

End Subroutine LerAij

Subroutine Metropolis()
    Use Variaveis
    Implicit None
    Integer :: J,i
    Real(8):: dE,r,esp
    Do j = 1,Ns
        i = j !int(rand()*Ns) + 1
        dE = -2.d0*S(i)*Bi(i)
        
        r = rand()
        !call random_number(r)
        Esp = exp(-beta*De)
        If (r.lt. esp) Then 

        If (beta > 1.d0 .and. beta < 4.d0 .and. dE>10) THEN
            Print*,i,de,esp,r
            Call config(beta)
            S(i)=-S(i)
            Call config(beta)
            S(i)=-S(i)
            Call flush()
            !pause
        End If


            Call Update(i,dE)
        End if
    End do         
End Subroutine Metropolis 


Subroutine MC()
    Use Variaveis
    Implicit None
    Integer :: i
    Real(8), Dimension(:), Allocatable :: Energia,Magx,Magy
    Real(8) :: Med, Var, Cv, med2,medx,medy
    character(70)::arq

    allocate(Energia(Nmc),Magx(Nmc),Magy(Nmc)) 

    Call Termaliza
    Do i = 1, Nmc
        Call Metropolis
        !Call Magnetizacao

        Energia(i) = Etot
        Magx(i)= Mx
        Magy(i)= My
    End Do 

    !Write(arq,"('energias_T=',F5.3,'.dat')")1.d0/Beta
    !Open(666,file=arq)
    !Do i=1,Nmc
    !    Write(666,*)i,energia(i)
    !End do
    !Close(666)


    Call Avevar(Nmc,Energia,med,med2,Var)
    deallocate(Energia)

    Cv = (beta**2/Real(Ns,8))*(med2 - med**2)

    Medx=sum(Magx)/dble(Ns)/dble(Nmc)
    Medy=sum(Magy)/dble(Ns)/dble(Nmc)

    Write(30,*) 1.d0/Beta, Med/dble(Ns), Cv, Medx, Medy
    Call Flush()

End Subroutine MC

Subroutine Avevar(n,Data,med,med2,var)
    Use Variaveis
    Implicit None
    integer, intent(in) :: n
    Real(8), intent(in) :: Data(n)    !! Também não muda valor
    Real(8) :: med, med2, var
    Real(8), dimension(:), allocatable :: Data2
    integer :: i
    Real(8) :: sp,ep

    med = 0.d0
    med = sum(data)/Real(N,8)
    allocate(Data2(n))

    Data2 = Data**2 ! Meu Data é a Energia 

    med2 = 0.d0
    med2 = sum(Data2)/Real(N,8)


    Var = 0.d0
    ep = 0.d0

    Do i = 1,n
        sp = data(i) - med
        ep = ep + sp
        Var = Var + Sp*Sp
    End Do
    
    Var = (Var - ep**2/Real(N,8))/Real(N-1,8)

    return

End Subroutine Avevar

Subroutine Temperatura()
    Use Variaveis
    Implicit None
    Integer :: i
    Integer :: nt = 15
    Real(8), dimension(:), allocatable :: T
    Real(8) :: Tf = 10.0d0, Ti = 50.d0, dT
    
    dT = 0.5

    Allocate(T(0:Nt))

    T(0)=Ti*2
        
    Do i = 1,Nt
        T(i) = T(i-1)*dt
    End Do

    Open(30,file = 'Emed.dat') 

    Do i = 1, Nt
        print*,T(i)
        Beta = 1.d0/T(i)
        !Print*, beta
        Call MC
        Call config(T(i))
    End Do   

    Close(30)

End Subroutine Temperatura

Subroutine Termaliza()
    Use Variaveis
    Implicit None
    Integer :: i

    Do i = 1,Nterma
        Call Metropolis
    End Do
    Return

End Subroutine Termaliza


Subroutine Inicia_BI()
    Use Variaveis
    Implicit None
    Integer :: I, J
    Allocate(Bi(Ns),S(Ns),Ei(Ns)) ! Spin Ising
    !S = 1
    
    Do i =1, Ns
        If (rand()<0.6) Then
            S(i)= 1
            Else
            S(i)= -1
        End If
    End Do
    
    Bi = 0.d0

    Do i = 1,Ns
        Do j = 1,Ns
            Bi(i) = Bi(i) + S(j)*Aij(i,j)
        End Do
        Ei(i) = S(i)*Bi(i) ! Energia em todos os spins  
    End Do
    Etot=sum(Ei)*0.5d0

    Mx=0.d0
    My=0.d0
    Mx=sum(ex*S)
    My=sum(ey*S)
    !Print*, 0.5d0*Sum(Ei)/Real(Ns,8)   ! Multiplicou por 1/2 porque conta duas vezes. (Energia Média por Spin)
End Subroutine Inicia_BI

Subroutine Update(i,dE)
    Use Variaveis
    Implicit None
    Integer, Intent(in) :: i !Não muda o valor de i.
    Real(8), Intent(in) :: dE! Não muda o valor de dE.
    Integer :: K, J
    Real(8) :: dBi

    Do k = 1,Ns
        dBi = -2.d0*S(i)*Aij(k,i)
        Bi(k) = Bi(k) + dBi
    End Do 
    S(i) = - S(i)
    Etot = Etot + dE

    Mx=Mx+2.d0*S(i)*Ex(i)
    My=My+2.d0*S(i)*Ey(i)

End Subroutine Update

Subroutine Config(t)
    Use Variaveis
    Implicit None
    Integer :: i
    Real(8), Dimension(:), allocatable :: Rx, Ry, Mx1, My1
    Real(8)::t
    Allocate(Rx(Ns),Ry(Ns),Mx1(Ns),My1(Ns))
    Open(10,file = 'resultados.csv')
    Read(10,*)
    Do i = 1,Ns
        Read(10,*)Rx(i),Ry(i),Mx1(i),My1(i)
    End Do
    Close(10)
    Open(4,file = 'Config.xyz')
    Write(4,*) Ns
    Write(4,*) Etot/Real(Ns,8),t
    Do i = 1,Ns
        Write(4,*)Rx(i),Ry(i),Mx1(i)*S(i),My1(i)*S(i), atan2(My1(i)*S(i),Mx1(i)*S(i)) ! ângulo do spin para colorir.
    End Do
    !Close(4)
End Subroutine Config

Program Main
    Use Variaveis
    Integer :: I
    Integer :: Seed = 387
    Call Srand(Seed)
    Call LerAij
    Call Inicia_BI
    Nterma = 10000
    Nmc = 100000
    !Call MC
    !Call Config
    Call Temperatura
End Program Main
