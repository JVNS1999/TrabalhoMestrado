!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                  Programa Cálculo do Campo Magnético                   !         
!                     Programador: João Vitor Nunes                      !                 
!                        Orientador: Lucas Mól                           !               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM field_calculation_and_plot
USE IFPORT
IMPLICIT NONE
REAL(KIND=8),PARAMETER:: delta=1.d-1, eps=1.d-2, length=7.d-1
REAL(8),DIMENSION(:),ALLOCATABLE:: rx, ry, psx, psy, sx, sy
INTEGER(4),DIMENSION(:),ALLOCATABLE:: flipind
REAL(8):: x, y, z, r, dx, dy, minx, maxx, miny, maxy
REAL(8):: nx, ny, r3, bx, by
INTEGER(4):: i, j, k, cont, nspins, npoints, npoints_real, nflip
LOGICAL(4):: exclog
CHARACTER(80):: cmd1
CALL RANDOM_SEED()

cmd1='wc -l < spins_pos_conf.dat > wcoutput.txt'
k = SYSTEM(cmd1)
OPEN(123,FILE='wcoutput.txt')
READ(123,*)nspins
CLOSE(123)
WRITE(*,'(A8,I6)')"nspins = ", nspins

ALLOCATE(psx(nspins),psy(nspins),sx(nspins),sy(nspins))

OPEN(1,FILE="spins_pos_conf.dat",STATUS='old')
DO i=1, nspins
        READ(1,*) k, psx(i), psy(i), sx(i), sy(i)
ENDDO
CLOSE(1) 

WRITE(*,*)"minx = ",MINVAL(psx)
WRITE(*,*)"maxx = ",MAXVAL(psx)
WRITE(*,*)"miny = ",MINVAL(psy)
WRITE(*,*)"maxy = ",MAXVAL(psy)
z = 2.d0
minx = MINVAL(psx)-z
maxx = MAXVAL(psx)+z
miny = MINVAL(psy)-z
maxy = MAXVAL(psy)+z

x = (maxx-minx)/delta
y = (maxy-miny)/delta
npoints = (IDINT(x)+1)*(IDINT(y)+1)
WRITE(*,*)"temporary npoints = ", npoints
ALLOCATE(rx(npoints),ry(npoints))
rx=-1.d4
ry=-1.d4
x = minx
y = miny
cont = 0
DO WHILE(x.LT.maxx)
        y = miny
        DO WHILE(y.LT.maxy)
                exclog=.FALSE.
                DO i=1, nspins
                        dx = x - psx(i)
                        dy = y - psy(i)
                        r = DSQRT(dx*dx+dy*dy)
                        IF(r.LT.eps)    THEN
                                exclog=.TRUE.
                        ENDIF
                ENDDO
                IF(exclog)      THEN
                        y=y+delta
                ELSE
                        cont=cont+1
                        rx(cont)=x
                        ry(cont)=y
                        y = y+delta
                ENDIF
        ENDDO
        x = x+delta
ENDDO
WRITE(*,*)cont
npoints_real=cont

OPEN(2,FILE="magfield.dat")
OPEN(3,FILE="rxry.dat")
OPEN(4,FILE="spins.dat")

nflip=5
ALLOCATE(flipind(nflip))
flipind(1)=308
flipind(2)=797
flipind(3)=409
flipind(4)=407
flipind(5)=398

DO i=1, npoints_real
        bx=0.d0
        by=0.d0
        DO j=1, nflip
                x=rx(i)-psx(flipind(j))
                y=ry(i)-psy(flipind(j))
                r=DSQRT(x*x+y*y)
                nx=x/r
                ny=y/r
                r3=1.d0/(r*r*r)
                bx=bx + (3.d0*nx*(sx(flipind(j))*nx+sy(flipind(j))*ny) - sx(flipind(j)))*r3
                by=by + (3.d0*ny*(sx(flipind(j))*nx+sy(flipind(j))*ny) - sy(flipind(j)))*r3
        ENDDO
        WRITE(2,'(2F14.6)')bx, by
        WRITE(3,'(2F14.6)')rx(i), ry(i)
ENDDO

CLOSE(4)
CLOSE(3)
CLOSE(2)

OPEN(5,FILE="sxsy_ini.dat")
OPEN(6,FILE="sxsy_end.dat")
OPEN(7,FILE="sxsy_flip_ini.dat")
OPEN(8,FILE="sxsy_flip_end.dat")

DO i=1,nspins
        z=DATAN2(sy(i),sx(i))
        x=psx(i)-0.5d0*length*DCOS(z)
        y=psy(i)-0.5d0*length*DSIN(z)
        WRITE(5,'(2F14.6)')x, y
        x=psx(i)+0.5d0*length*DCOS(z)
        y=psy(i)+0.5d0*length*DSIN(z)
        WRITE(6,'(2F14.6)')x,y
ENDDO
DO i=1,nflip
        z=DATAN2(sy(flipind(i)),sx(flipind(i)))
        x=psx(flipind(i))-0.5d0*length*DCOS(z)
        y=psy(flipind(i))-0.5d0*length*DSIN(z)
        WRITE(7,'(2F14.6)')x, y
        x=psx(flipind(i))+0.5d0*length*DCOS(z)
        y=psy(flipind(i))+0.5d0*length*DSIN(z)
        WRITE(8,'(2F14.6)')x,y
ENDDO



CLOSE(8)
CLOSE(7)
CLOSE(6)
CLOSE(5)

OPEN(555,FILE="num_lines.dat")
WRITE(555,*)npoints_real
WRITE(555,*)nspins
WRITE(555,*)nflip
CLOSE(555)

DEALLOCATE(rx,ry,psx,psy,sx,sy,flipind)

STOP
END PROGRAM field_calculation_and_plot
