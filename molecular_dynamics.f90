program molecular_dynamics
        !init, force and integrate son subrutinas 

        implicit none

        save            ! le dice a las subrutinas declaradas que trabajen al mismo tiempo

        integer :: npart, nx, ny, nz
        !real(8) significa que las variables tendrán el doble de precision (6 dig ---> 12 dig)
        real(8) :: delt, t, tmax, Lx, Ly, Lz, temp, rcut
        real(8) :: ecin, epot, etot
        real(8), allocatable, dimension(:,:)  :: x, v, xm, f 

        call read_data
        
        !write(*,*)delt
        
        call init
        
        t = 0.0

        do while (t <= tmax)

                call force !(f,n)
                call integrate !(f,n)

                t = t + delt

!                call sample

        !        call pbc        periodic_boundory_condition 

                call jmol_visualize

        end do

        stop

        contains        !se guardan los procedimientos auxiliares en este bloque

        ! *********************************************************************!
        !       AQUI SE LEEN LOS DATOS DE LA SIMULACION                        !
        ! *********************************************************************!

        subroutine read_data

                implicit none

                write(*,*)'Leyendo los parametros de simulacion'

                open(100,file="input_data.inp")
                        read(100,*)
                        read(100,*)nx, ny, nz   !numero de atomos en cada direccion
                        read(100,*)
                        read(100,*)Lx, Ly, Lz   !longitud de los lados de  la caja de simulacion
                        read(100,*)
                        read(100,*)temp
                        read(100,*)
                        read(100,*)tmax, delt
                        read(100,*)
                        read(100,*)rcut
                close(100)

                write(*,*)'Los parametros han sido leidos'

                npart = nx*ny*nz

                write(*,*)"Alojando arreglos de posicion y velocidad"
                write(*,*)"para",npart," atomos"

                allocate(x(1:npart,1:3), v(1:npart,1:3), xm(1:npart,1:3), f(1:npart,1:3))

        end subroutine read_data 

       !**********************************************************************!
       !                INICIALIZACION                                        !
       !**********************************************************************!

       subroutine init

               implicit none

               integer :: i, j, k, ipart        !ipart es la i-esima particula
               real(8) :: delx, dely, delz, sumv2, fs
               real(8), dimension(1:3) :: sumv

               delx = Lx/real(nx)
               dely = Ly/real(ny)
               delz = Lz/real(nz)

               ipart = 1

               write(1000,*)npart
               write(1000,*)"mi simulacion"

               do k = 1, nz
                        do j = 1, ny
                                do i = 1, nx
                                        
                                        x(ipart,1) = 0.5*delx + real(i-1)*delx
                                        x(ipart,2) = 0.5*dely + real(j-1)*dely
                                        x(ipart,3) = 0.5*delz + real(k-1)*delz
                                        
                                        write(1000,*)"Au",x(ipart,1),x(ipart,2),x(ipart,3)

                                        ipart = ipart + 1

                                end do
                        end do
               end do
               
               ! Inicialización de velocidades

               sumv(1) = 0.0
               sumv(2) = 0.0
               sumv(3) = 0.0

               sumv2 = 0.0

               do i = 1, npart
               
                        ! velocidades iniciales con distribucion uniforme

                        v(i,1) = -0.5 + rand(0) 
                        v(i,2) = -0.5 + rand(0)
                        v(i,3) = -0.5 + rand(0)

                        !velocidad del centro de masa
                        sumv(1) = sumv(1) + v(i,1)
                        sumv(2) = sumv(2) + v(i,2)
                        sumv(3) = sumv(3) + v(i,3)

                        !energia cinetica inecial

                        sumv2 = sumv2 + v(i,1)**2 + v(i,2)**2 + v(i,3)**2 

               end do

               sumv(1) = sumv(1)/real(npart)
               sumv(2) = sumv(2)/real(npart)
               sumv(3) = sumv(3)/real(npart)

               fs = sqrt(3.0*real(npart)*temp/sumv2)

               !restitucion de la velocidad del centro de masa y reescalamiento

               do i = 1, npart
                        
                        v(i,1) = fs * (v(i,1) - sumv(1))
                        v(i,2) = fs * (v(i,2) - sumv(2))
                        v(i,3) = fs * (v(i,3) - sumv(3))

                        !primeras posiciones previas como MRU

                        xm(i,1) = x(i,1) - v(i,1) + delt
                        xm(i,2) = x(i,2) - v(i,2) + delt
                        xm(i,3) = x(i,3) - v(i,3) + delt

               end do

       end subroutine init

       !**********************************************************************!
       !                   CALCULO DE LAS FUERZAS                             !     
       !**********************************************************************!

       subroutine force

               implicit none

               integer :: i, j
               real(8) :: xrel, yrel, zrel, dist2, rcut2, fx, fy, fz

               ! Inicialización de las fuerzas 

               !f es una matrix, pero con la orden que le asigna el valor 0.0,
               !estamos forzando a que Fortran haga cero todas las entradas de f
               !Se trata de una operacion compacta y funciona con otro valor
               !distinto de cero, incluso podemos compactar de esta manera la 
               !multiplicacion de matrices.
               f = 0.0 

               rcut2 = rcut**2

               do i = 1, npart - 1
                        do j = i+1, npart
                                !recuerda que la primera entrada de x hace referencia a la componente mientras que
                                !la segunda entrada hace referencia a los ejes x,y,z
                                xrel = x(i,1) - x(j,1)
                                yrel = x(i,2) - x(j,2)
                                zrel = x(i,3) - x(j,3)

                                !Transformacion a la imagen mas cercana

                                xrel = xrel - Lx * nint(xrel/Lx)
                                yrel = yrel - Ly * nint(yrel/Ly)
                                zrel = zrel - Lx * nint(zrel/Lz)

                                ! Distancia al cuadrado

                                dist2 = xrel**2 + yrel**2 + zrel**2

                                !Consideracion del radio de corte
                                
                                !rcut sera el radio de corte
                                if (dist2 < rcut2)then
                                        
                                        fx = (48.0/dist2)*((1/dist2**6)- 0.5*(1/dist2**3))*xrel 
                                        fy = (48.0/dist2)*((1/dist2**6)- 0.5*(1/dist2**3))*yrel 
                                        fz = (48.0/dist2)*((1/dist2**6)- 0.5*(1/dist2**3))*zrel 

                                        f(i,1) = f(i,1) + fx
                                        f(i,2) = f(i,2) + fy
                                        f(i,3) = f(i,3) + fz
 
                                        f(j,1) = f(j,1) - fx
                                        f(j,2) = f(j,2) - fy
                                        f(j,3) = f(j,3) - fz
                                        
                                      !  write(*,*)fx, fy, fz
                                end if
                        end do
               end do

       end subroutine force

       !**********************************************************************!
       !     INTEGRACION DE LAS ECUACIONES DE MOVIMIENTO                      !
       !**********************************************************************!

       subroutine integrate

                implicit none

                integer :: i
                real(8) :: sumv2
                real, dimension(1:3) :: vcm, s

                sumv2 = 0.0 
                vcm   = 0.0

                do i = 1, npart
                ! Aplicación de la fórmula de Verlet
                s(1) = 2.0*x(i,1) - xm(i,1) + f(i,1) * delt**2
                s(2) = 2.0*x(i,2) - xm(i,2) + f(i,2) * delt**2
                s(3) = 2.0*x(i,3) - xm(i,3) + f(i,3) * delt**2

                v(i,1) =0.5*(s(1) -xm(i,1))/delt
                v(i,2) =0.5*(s(2) -xm(i,2))/delt
                v(i,3) =0.5*(s(3) -xm(i,3))/delt

                !contribucion a la velocidad del centro de masa

                vcm(1) = vcm(1) + v(i,1)
                vcm(2) = vcm(2) + v(i,2)
                vcm(3) = vcm(3) + v(i,3)

                !contribucion a la energia cinetica

                sumv2 = sumv2 + v(i,1)**2 + v(i,2)**2 + v(i,3)**2

                !actualizacion de las posiciones

                !posicion previa

                xm(i,1) = x(i,1)
                xm(i,2) = x(i,2)
                xm(i,3) = x(i,3) 

                !posicion actual
 
                x(i,1) = s(1)
                x(i,2) = s(2)
                x(i,3) = s(3)

                end do

       vcm = vcm / real(npart)
       ecin = sumv2/(2.0*real(npart))
       temp = sumv2/(3.0*real(npart))

       write(2001,*)t, vcm(1),  vcm(2), vcm(3)
       write(2002,*)t, ecin
       write(2003,*)t, temp

       end subroutine integrate 

       !**********************************************************************!
       !         ESCRITURA PARA VISUALIZACION EN JMOL                         !
       !**********************************************************************!
       
       subroutine jmol_visualize

               implicit none

               integer :: i 
               character :: atname, mssg

               atname = "Au"
               mssg = "My simulation"

               write(1000,*)npart
               write(1000,*)mssg

               do i = 1, npart

                        write(1000,*)"Au", x(i,1), x(i,2), x(i,3)
         
               end do
               
       end subroutine jmol_visualize 

       
       !subroutine pbc
                !IMPLEMENTACION DE CONDICIONES DE FRONTERA PERIODICAS

        !        do i = 1, npart
         !               x(i,1) = x(i,1) - Lx * floor(x(i,1)/Lx)
          !              x(i,2) = x(i,2) - Lx * floor(x(i,2)/Lx)
           !             x(i,3) = x(i,3) - Lx * floor(x(i,3)/Lx)

                        
            !    end do
      ! end subroutine pbc
end program molecular_dynamics
