Module Variaveis
    Implicit None
    Integer :: Ns, Nterma
    Integer :: Nmc
    Integer, Dimension(:), allocatable :: S
    Real(8), Dimension(:,:), allocatable :: Aij
    Real(8), Dimension(:),allocatable :: Bi, Ei
    Real(8) :: beta , Etot
End Module Variaveis

Subroutine LerAij()
    Use Variaveis
    Implicit None
    Integer :: I , J
    Open(20,file = 'Aij.dat', status = 'Old')
    Read(20,*)Ns
    Allocate(Aij(Ns,Ns))
    Do i = 1,Ns
        Do j = 1,Ns
            Read(20,*) Aij(i,j)
        End Do
    End Do
    Close(20)     
End Subroutine LerAij

Subroutine Metropolis()
    Use Variaveis
    Implicit None
    Integer :: J,i
    Real(8):: dE
    Do j = 1,Ns
        i = int(rand()*Ns) + 1
        dE = -2.d0*S(i)*Bi(i)
        If (rand().lt. exp(-beta*dE)) then 
            Call Update(i,dE)
        End if
    End do         
End Subroutine Metropolis 


Subroutine MC()
    Use Variaveis
    Implicit None
    Integer :: i
    Real(8), Dimension(:), Allocatable :: Energia 
    Real(8) :: Med, Var, Cv, med2
    allocate(Energia(Nmc)) 

    Call Termaliza
    Do i = 1, Nmc
        Call Metropolis
        Energia(i) = Etot
    End do 

    Call Magnet(Nmc,Energia,med,med2,Var)
    deallocate(Energia)

    Cv = (beta**2/Real(Ns,8))*(med2 - med**2)
    Write(30,*) 1.d0/Beta, Med, Var, Cv
    Call Flush()

End Subroutine MC

Subroutine Magnet(n,Data,med,med2,var)
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

End Subroutine Magnet

Subroutine Temperatura()
    Use Variaveis
    Implicit None
    Integer :: i
    Integer :: nt = 50
    Real(8), dimension(:), allocatable :: T
    Real(8) :: Tf = 0.01d0, Ti = 2.d0, dT
    
    dT = (Tf - Ti)/Real(nt,8)

    Allocate(T(nt))

    Do i = 1,nt
        T(i) =  Ti + (i-1)*dT
    End Do

    Open(30,file = 'Emed.dat') 

    Do i = 1, Nt
        Beta = 1.d0/T(i)
        Print*, beta
        Call MC
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
    S = 1
    Bi = 0.d0
    Do i = 1,Ns
        Do j = 1,Ns
            Bi(i) = Bi(i) + S(j)*Aij(i,j)
        End Do
        Ei(i) = S(i)*Bi(i) ! Energia em todos os spins  
    End Do
    !Print*, 0.5d0*Sum(Ei)/Real(Ns,8)   ! Multiplicou por 1/2 porque conta duas vezes. (Energia Média por Spin)
End Subroutine Inicia_BI

Subroutine Update(i,dE)
    Use Variaveis
    Implicit None
    Integer, Intent(in) :: i !Não muda o valor de i.
    Real(8), Intent(in) :: dE ! Não muda o valor de dE.
    Integer :: K, J
    Real(8) :: dBi
    Do k = 1,Ns
        dBi = -2.d0*S(i)*Aij(k,i)
        Bi(k) = Bi(k) + dBi
    End Do 
    S(i) = - S(i)
    Etot = Etot + dE  
End Subroutine Update

Subroutine Config()
    Use Variaveis
    Implicit None
    Integer :: i
    Real(8), Dimension(:), allocatable :: Rx, Ry, Mx, My
    Allocate(Rx(Ns),Ry(Ns),Mx(Ns),My(Ns))
    Open(10,file = 'resultados.csv')
    Read(10,*)
    Do i = 1,Ns
        Read(10,*)Rx(i),Ry(i),Mx(i),My(i)
    End Do
    Close(10)
    Open(4,file = 'Config.xyz')
    Write(4,*) Ns
    Write(4,*) Etot/Real(Ns,8)
    Do i = 1,Ns
        Write(4,*)Rx(i),Ry(i),Mx(i)*S(i),My(i)*S(i), atan2(My(i)*S(i),Mx(i)*S(i)) ! ãngulo do spin para colorir.
    End Do
    !Close(4)
End Subroutine Config

Program Main
    Use Variaveis
    Integer :: I
    Integer :: Seed = 1234
    Call Srand(Seed)
    Call LerAij
    Call Inicia_BI
    Nterma = 10000
    Nmc = 1000000
    !Call MC
    !Call Config
    Call Temperatura
End Program Main
