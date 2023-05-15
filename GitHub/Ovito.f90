Program MonteCarlo
 Implicit None 
 
  !Declaração de Variáveis
  Real(8), Dimension(:),allocatable :: x_spin, y_spin, xs, ys ! Dimensão e Alocada na Memória !!
  Integer, Dimension(:),allocatable :: S
  Integer:: Ns = 275, Seed = 1234
  Integer:: I, J, K
  Character(10):: Dummy
  
  Open(10, file = 'resultados.csv', status = 'Old')
  Read(10,*)
  Allocate(x_spin(Ns), y_spin(Ns), xs(Ns), ys(Ns),S(Ns))
  S = 1 ! Todos os valores de S são unitários.
    Do i = 1,Ns
     Read(10, *) x_spin(i), y_spin(i), xs(i), ys(i)
    End Do
  Close(10)

  Call SRand(Seed) ! Semente de gerador de número aleatório.
  ! Print*, x_spin, y_spin, xs, ys

   Open(20, file ='Rede.xyz')
    Do J = 1,10
   Write(20, *) Ns 
   Write(20, *) "C"
    Do i = 1,Ns
      Write(20, *) x_spin(i), y_spin(i), S(i)*xs(i), S(i)*ys(i), atan2(S(i)*ys(i), S(i)*xs(i))
    End Do
    Do i = 1, Ns
     K = int(Ns*rand()) + 1 ! Multiplico o número entre 0 e 1 - Array, rand é um número real multiplicando por Ns eu escalo. A soma vem dos limites do Array.
     If(Rand() > 0.5) then
     S(k) = -S(k)
     End If
     End Do
    End Do
    Close(20)

End Program MonteCarlo