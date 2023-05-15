Module Rede
    Integer, Parameter :: Nvertice = 381 ! Número de Vértices
    Integer:: Nspin !Número de Spin
    Real(8),Dimension(:), Allocatable :: Vx,Vy, x, y, mx, my
End Module    


Subroutine Ler_Vertices
    Use Rede
    Implicit None
    Integer :: I
    Real(8):: Z
    
    Open(10, file = 'tiling.csv')
    Open(12, file = 'Vertice.xyz')
    Write(12,*) Nvertice
    Write(12,*) ''  
    Allocate(Vx(Nvertice),Vy(Nvertice))
    Do i = 1, Nvertice
        Read(10,*) Vx(i),Vy(i), Z
        Write(12,*) Vx(i),Vy(i)
    End Do
    Close(10)
    Close(12)
End Subroutine

Subroutine Criar_Rede
    Use Rede
    Implicit None
    Integer :: i, j, k, ii
    Real(8) :: dx, dy, dist


    Open(10, file = 'rede.dat')
    K = 0 !Contador
    !VSV = 0 !Vértice-Spin-Vértice
    !Allocate(VSV(Nvertice,Nspin))
    Do i = 1,Nvertice - 1
        Do j = i + 1, Nvertice
            dx = Vx(i) - Vx(j)
            dy = Vy(i) - Vy(j)
            dist = sqrt(dx**2 + dy**2)
            If(dist < 1.01d0 .and. dist > 0.99d0) Then
                K = K + 1
                Write(10,*) Vx(i) - 0.5d0*dx, Vy(i) - 0.5d0*dy, dx/dist, dy/dist
                !x(K) = Vx(i) - 0.5d0*dx
                !y(K) = Vy(i) - 0.5d0*dy
                !mx(K) = dx/dist
                !my(K) = dy/dist
                !VSV(i,k) = j
                !VSV(j,k) = i
            End If
        End Do
    End Do
    Close(10)
    Nspin = K
    Allocate(x(Nspin), y(Nspin), mx(Nspin), my(Nspin))
    Open(10, file = 'rede.dat')

    Do i = 1,Nspin
        Read(10,*) x(i),y(i), mx(i), my(i)
    End Do     
    Close(10)

    Open(11, file='config.xyz')
    Write(11,*) Nspin
    Write(11,*) ''
    Do i = 1,Nspin
        Write(11,*) x(i),y(i), mx(i), my(i)
    End Do
    Close(11)

    Return
End Subroutine Criar_Rede

Subroutine Vizinhos
    Use Rede
    Implicit None
    Integer :: i, j, k 
    Real(8) :: dx,dy,dist,Jij

    Open(10, file = 'VizVertice.dat')
    Open(11, file = 'Nviz.dat')
    Write(11,*) 0 !Recebe um a mais!
    K = 0
    Do i = 1,Nvertice
        Do j = 1,Nspin
            dx = Vx(i) - x(j) ! Vértice - Spin
            dy = Vy(i) - y(j) ! Vértice - Spin
            dist = sqrt(dx**2 + dy**2)
            If (dist < 0.6d0 .and. dist> 0.4d0) Then
                K = K + 1
                Jij = dx*mx(j) + dy*my(j)
                Write(10,*) i, j, Int(Sign(1.d0, Jij))
            End If
        End Do
        Write(11,*)K
    End Do
    Close(10)
    Close(11)
    Return    
End Subroutine Vizinhos
!! Apontando para o vértice é positivo, saindo é negativo.

Program Main
    Call Ler_Vertices
    Call Criar_Rede
    Call Vizinhos
End Program Main