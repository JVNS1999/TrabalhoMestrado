Module Rede
    Implicit None
    Integer :: Ns = 275
    Real(8), Dimension(:), allocatable :: Rx, Ry, Mx, My
End Module Rede

Subroutine Inicial()
    Use Rede
    Implicit None
    Integer :: I
    

    Allocate(Rx(Ns),Ry(Ns),Mx(Ns),My(Ns))
    Open(10,file = 'resultados.csv')
    Read(10,*)
    Do i = 1,Ns
        Read(10,*)Rx(i),Ry(i),Mx(i),My(i)
    End Do
    Close(10)
    Call CalcCampo  
End Subroutine Inicial 

Subroutine CalcCampo()
    Use Rede
    Implicit None
    Integer :: I , J
    Real(8) :: dx, dy, dist, E1, E2, Eij
    Real(8), Dimension(:,:),allocatable :: Aij
    Allocate(Aij(Ns,Ns - 1))
    Do I = 1,Ns
        Do J = 1,Ns
            If (I .ne. J) then
                dx = Rx(i) - Rx(j)
                dy = Ry(i) - Ry(j)
                dist = sqrt(dx**2 + dy**2)
                dx = dx/dist 
                dy = dy/dist
                E1 = mx(i)*mx(j) + my(i)*my(j)
                E2 = -3.d0*(mx(i)*dx + my(i)*dy)*(mx(j)*dx + my(j)*dy)
                Eij = (E1 + E2)/dist**3
                Aij(i,j) = Eij ! Salvando numa matriz os valores de interação.
            End If
        End Do        
    End Do
    Open(20, file = 'Aij.dat')
    Write(20,*) Ns
    Do i = 1,Ns
        Do j = 1, Ns - 1 ! Não tenho interação com o mesmo spin por isso o -1.
            If (i .ne. j) then 
                Write(20,*) Aij(i,j)
            End If    
        End Do
    End Do
    Close(20)        
        
End Subroutine CalcCampo    
Program Main
    Call Inicial
End Program Main