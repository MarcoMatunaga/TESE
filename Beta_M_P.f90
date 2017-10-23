module vars
implicit none
!
!
real(8),dimension(:,:),allocatable               :: meshx, meshy
real(8),dimension(:,:),allocatable               :: residue1,residue2,residue3,residue4
real(8),dimension(:,:),allocatable               :: E_flux1,E_flux2,E_flux3,E_flux4 
real(8),dimension(:,:),allocatable               :: F_flux1,F_flux2,F_flux3,F_flux4 
real(8),dimension(:,:),allocatable               :: E_Harten1,E_Harten2,E_Harten3,E_Harten4 
real(8),dimension(:,:),allocatable               :: F_Harten1,F_Harten2,F_Harten3,F_Harten4 
real(8),dimension(:,:),allocatable               :: u_vetor, v_vetor, p_vetor, c_vetor, H_vetor
real(8),dimension(:,:),allocatable               :: u_roe, v_roe, p_roe, c_roe, H_roe
real(8),dimension(:),allocatable                 :: R1,R2,R3,R4  
real(8),dimension(:),allocatable                 :: x,y  
real(8)                                          :: delta_rho, delta_e, delta_rhou,delta_rhov
real(8)                                          :: deltat,Deltax, Deltay, Height, Length 
real(8)                                          :: aux3, aux4, C1, C2, alfa1, alfa2,alfa3,alfa4
real(8)                                          :: p,u,v,rho,gama,residue, tempo_total
real(8)                                          :: eps1,eps2,eps3,eps4, a,A1
real(8)                                          :: eigen1,eigen2,eigen3,eigen4
real(8)                                          :: arg1,arg2,arg3,arg4
real(8)                                          :: max_residue, res, valor, CFL
real(8)                                          :: shock_angle, mach_inlet, turn_angle, mach_inlet_n
real(8)                                          :: mach_2, mach_2_n, trans, mach_core, p_percent
integer(4)                                       :: i,j, points_x, points_y, i_n
integer(4)                                       :: iter, dimen, max_iter, backup, arq, tec_sol
real(8),dimension(:,:),allocatable               :: Q1,Q2,Q3,Q4
real(8),dimension(4)                             :: psi
!
!
    contains 
            !
        real(kind=8) function q_z(z,eps)
            !
            implicit none
            !
            real(kind=8) :: z, eps
            !
            q_z = dabs(z)
            !
            if (dabs(z) < eps) then
                q_z = 0.5d0*( (z**2.0d0/eps) + eps)
            end if
            !
        end function q_z
!
!
end module vars
!
!**********Programa_projeto4**********
!
program proj4
use vars
implicit none
!
! Dados do problema
!
namelist /PAR_Flow/ mach_inlet,p_percent,shock_angle,gama
namelist /PAR_Time/ max_iter,tempo_total
namelist /PAR_Mesh/ points_x,points_y,Height,Length
!
!
!
open(11,file='inputs')
read(11,PAR_Flow)
read(11,PAR_Time)
read(11,PAR_Mesh)
close(11)
!
!
!
!mach_core    = sqrt(2.0d0*(p_percent)/gama)
A1 = ( (p_percent/100.0d0)**((gama -1.0d0)/gama))*(1.0d0 + ((gama -1.0d0)/2.0d0)*mach_inlet**2.0d0 )
mach_core    = sqrt( (A1 - 1.0d0)*2.0d0/(gama-1.0d0) )
shock_angle  = shock_angle*(3.1415d0/180.0d0)
!
! calculo para determinar as condicoes pos choque
!
turn_angle   = datan( (2.0d0/tan(shock_angle)) * ( (mach_inlet*sin(shock_angle))**2.0d0 - 1.0d0 )/( (mach_inlet**2.0d0)*(gama + cos(2.0d0*shock_angle)) + 2.0d0)  )
!turn_angle   = 10.94d0*(3.1415d0/180.0d0)
mach_inlet_n = mach_inlet*sin(shock_angle) 
mach_2_n     = sqrt((2.0d0 + (gama - 1.0d0)*mach_inlet_n**2.0d0)/( 2.0d0*gama*mach_inlet_n**2.0d0 - (gama - 1.0d0) ) )
mach_2       =  mach_2_n/sin(shock_angle - turn_angle)
!
! dimensoes vetores de fluxo e de variaveis
!
dimen = 4
!
! marcha no tempo
! 
deltat = tempo_total/max_iter
iter = 0
tec_sol = 1
arq = 0
!
! Criterio de convergencia
!
max_residue = 1e-5
residue = 10
!
! Delta x e y, sao respectivamente os espacamentos de malha
! em x e y 
!
Deltax = Length/DBLE(points_x)
Deltay = Height/DBLE(points_y)
trans = 0.2d0/Deltay
i_n = int(trans)
!
!
allocate(E_flux1(points_x,points_y),E_flux2(points_x,points_y))
allocate(E_flux3(points_x,points_y),E_flux4(points_x,points_y))
allocate(E_Harten1(points_x,points_y),E_Harten2(points_x,points_y))
allocate(E_Harten3(points_x,points_y),E_Harten4(points_x,points_y))
allocate(F_flux1(points_x,points_y),F_flux2(points_x,points_y))
allocate(F_flux3(points_x,points_y),F_flux4(points_x,points_y))
allocate(F_Harten1(points_x,points_y),F_Harten2(points_x,points_y))
allocate(F_Harten3(points_x,points_y),F_Harten4(points_x,points_y))
allocate(u_vetor(points_x,points_y),v_vetor(points_x,points_y))
allocate(p_vetor(points_x,points_y),c_vetor(points_x,points_y))
allocate(u_roe(points_x,points_y),v_roe(points_x,points_y))
allocate(p_roe(points_x,points_y),c_roe(points_x,points_y))
allocate(H_vetor(points_x,points_y))
allocate(H_roe(points_x,points_y))
allocate(R1(dimen),R2(dimen))
allocate(R3(dimen),R4(dimen))
allocate(residue1(points_x,points_y),residue2(points_x,points_y))
allocate(residue3(points_x,points_y),residue4(points_x,points_y))
allocate(meshx(points_x,points_y), meshy(points_x,points_y))
allocate(x(points_x), y(points_y))
allocate(Q1(points_x,points_y), Q2(points_x,points_y))
allocate(Q3(points_x,points_y), Q4(points_x,points_y))
!
!
call mesh
call initial_cond
call boundary_conditions
!
!
!    .and. 
do while(iter < max_iter)
    call fluxes
    call fluxHarten_x
    call fluxHarten_y
    call convergence
    call Harten
    !call Harten_implicit
    !call Harten_implicit2
    call boundary_conditions
    iter = iter + 1
    backup = tec_sol*max_iter/20
    if (backup == iter .and. tec_sol <= 20)then
        call transient
        tec_sol = tec_sol + 1
    end if
!    open(1,file='residue.txt')
!    open(2,file='CFL.txt')
!    write(1,*) iter, log10(abs(residue))
!    write(2,*) iter, CFL
    !print *, iter
enddo
!    close(1)
!    close(2)
!
!
!
call output
!
!
deallocate(E_flux1,E_flux2)
deallocate(E_flux3,E_flux4)
deallocate(E_Harten1,E_Harten2)
deallocate(E_Harten3,E_Harten4)
deallocate(F_flux1,F_flux2)
deallocate(F_flux3,F_flux4)
deallocate(F_Harten1,F_Harten2)
deallocate(F_Harten3,F_Harten4)
deallocate(u_vetor,v_vetor)
deallocate(p_vetor,c_vetor)
deallocate(u_roe,v_roe)
deallocate(p_roe,c_roe)
deallocate(H_vetor)
deallocate(H_roe)
deallocate(R1,R2)
deallocate(R3,R4)
deallocate(residue1,residue2)
deallocate(residue3,residue4)
deallocate(x, y)
deallocate(meshx,meshy,Q1,Q2,Q3,Q4)
!
!
end program proj4
!
!**********mesh**********
!
subroutine mesh
use vars
implicit none
!
!
!
    do i = 1, points_x
        x(i) = i*Deltax
    end do

    do j = 1, points_y
        y(j) = j*Deltay
    end do

    do i = 1, points_x
        do j = 1, points_y
            meshx(i,j) = x(i)
            meshy(i,j) = y(j)
        end do
    end do
!
!
end subroutine mesh
!
!**********initial_conditions**********
!
subroutine initial_cond
use vars
implicit none
!
! Mach number 2.9/p 0.72atm/T 400K/
!
        p   = 1.0d0/gama
        rho = 1.0d0
        u   = mach_inlet
        v   = 0.0d0
        !
        !
    do i = 1, points_x
        do j =1, points_y/2
            Q1(i,j) =  rho
            Q2(i,j) =  rho*u
            Q3(i,j) =  rho*v
            Q4(i,j) =  p/(gama-1.0d0) + 0.50d0*rho*(u**2.0d0 + v**2.0d0)
        enddo
    enddo
    !
    !
        u   = mach_2
        v   = 0.0d0
        rho =  (gama + 1.0d0)*mach_inlet_n**2.0d0/(2.0d0 + (gama - 1.0d0)*mach_inlet_n**2.0d0)
        p   =  (1.0d0 + (2.0d0*gama/( gama + 1.0d0 ))*(mach_inlet_n**2.0d0 - 1.0d0) ) / gama
        !
        !
    do i = 1, points_x
        do j = points_y/2 + 1, points_y
            Q1(i,j) =  rho
            Q2(i,j) =  rho*u
            Q3(i,j) =  rho*v
            Q4(i,j) =  p/(gama-1.0d0) + 0.50d0*rho*(u**2.0d0 + v**2.0d0)
        enddo
    enddo
!
!
!
end subroutine initial_cond
!!
!!
subroutine boundary_conditions
use vars
implicit none
!
! entrada
!
    i = 1
    !
    ! Região 1
    !
        u   = mach_core
        v   = 0.0d0
        rho = 1.0d0
        p   = 1.0d0/gama
		!
		!
        do j = 1, i_n
            Q1(i,j) =  rho
            Q2(i,j) =  rho*u
            Q3(i,j) =  rho*v
            Q4(i,j) =  p/(gama-1.0d0) + 0.50d0*rho*(u**2.0d0 + v**2.0d0)
        enddo
    !
    ! Região 2
    !
        u   = mach_inlet
        v   = 0.0d0
        rho = 1.0d0
        p   = 1.0d0/gama
        !
        !
        do j = i_n+1, points_y/2
            Q1(i,j) =  rho
            Q2(i,j) =  rho*u
            Q3(i,j) =  rho*v
            Q4(i,j) =  p/(gama-1.0d0) + 0.50d0*rho*(u**2.0d0 + v**2.0d0)
        enddo
    !
    ! região 3
    !
        u   = mach_2
        v   = 0.0d0
        rho =  (gama + 1.0d0)*mach_inlet_n**2.0d0/(2.0d0 + (gama - 1.0d0)*mach_inlet_n**2.0d0)
        p   =  (1.0d0 + (2.0d0*gama/( gama + 1.0d0 ))*(mach_inlet_n**2.0d0 - 1.0d0) ) / gama
        !
        !
        do j = points_y/2 + 1, points_y
            Q1(i,j) =  rho
            Q2(i,j) =  rho*u
            Q3(i,j) =  rho*v
            Q4(i,j) =  p/(gama-1.0d0) + 0.50d0*rho*(u**2.0d0 + v**2.0d0)
        enddo
!
! froonteira inferior --> parede 
!
    do i = 1, points_x

        j = 1

        v = 0.0d0
        u = Q2(i,j+1)/Q1(i,j+1)
        rho = Q1(i,j+1)
        p = (gama-1.0d0) * (Q4(i,j+1) - 0.5d0*(Q2(i,j+1)**2.0d0 + Q3(i,j+1)**2.0d0)/rho)

        Q1(i,j) = rho
        Q2(i,j) = rho*u
        Q3(i,j) = rho*v
        Q4(i,j) = (p/(gama-1.0d0)) + 0.5d0*rho*(u**2.0d0 + v**2.0d0)

    end do
!
! saida 
!
    do j = 1, points_y

        i = points_x

        Q1(i,j) = Q1(i-1,j)
        Q2(i,j) = Q2(i-1,j)
        Q3(i,j) = Q3(i-1,j)
        Q4(i,j) = Q4(i-1,j)

    end do
!
! Fronteira superior --> parede
!       
    do i = 1, points_x
        
        j = points_y
        
        v = 0.0d0
        u = Q2(i,j-1)/Q1(i,j-1)
        rho = Q1(i,j-1)
        p = (gama-1.0d0) * (Q4(i,j-1) - 0.5d0*(Q2(i,j-1)**2.0d0 + Q3(i,j-1)**2.0d0)/rho)

        Q1(i,j) = rho
        Q2(i,j) = rho*u
        Q3(i,j) = rho*v
        Q4(i,j) = (p/(gama-1.0d0)) + 0.5d0*rho*(u**2.0d0 + v**2.0d0)

    end do

end subroutine
!
!**********fluxes**********
!
subroutine fluxes
use vars
implicit none
!
!
    do i = 1,points_x
        do j = 1,points_y
            !
            !
            u = Q2(i,j)/Q1(i,j)
            v = Q3(i,j)/Q1(i,j)
            p = (gama - 1.0d0)*(Q4(i,j) - 0.5d0*Q1(i,j)*(u**2.0d0 + v**2.0d0))
            !
            !
            E_flux1(i,j) = Q2(i,j)
            E_flux2(i,j) = Q1(i,j)*u**2.0d0 + p
            E_flux3(i,j) = Q1(i,j)*u*v
            E_flux4(i,j) = (Q4(i,j) + p)*u
            !
            !
        enddo
    enddo
    !
    !
    do j = 1,points_y
        do i = 1,points_x
            !
            !
            u = Q2(i,j)/Q1(i,j)
            v = Q3(i,j)/Q1(i,j)
            p = (gama - 1.0d0)*(Q4(i,j) - 0.5d0*Q1(i,j)*(u**2.0d0 + v**2.0d0))
            !
            !
            F_flux1(i,j) = Q3(i,j)
            F_flux2(i,j) = Q1(i,j)*u*v
            F_flux3(i,j) = Q1(i,j)*v**2.0d0 + p
            F_flux4(i,j) = (Q4(i,j) + p)*v
            !
            !
        enddo
    enddo
!
!
end subroutine fluxes
!
!**********fluxos de Harten na direcao x**********
!
subroutine fluxHarten_x
use vars
implicit none
!
!
real(8), dimension(:,:), allocatable            :: g1_x,g2_x,g3_x,g4_x
real(8), dimension(:,:), allocatable            :: gtil1_x,gtil2_x,gtil3_x,gtil4_x
real(8), dimension(:,:), allocatable            :: s1,s2,s3,s4
real(8), dimension(:,:), allocatable            :: gama_harten1_x,gama_harten2_x,gama_harten3_x,gama_harten4_x
real(8), dimension(:,:), allocatable            :: ni1_x,ni2_x,ni3_x,ni4_x 
!
! eps
!
    eps1 = 0.01d0
    eps2 = 0.0d0
    eps3 = 0.0d0
    eps4 = 0.0d0
!
!
allocate(g1_x(points_x,points_y),g2_x(points_x,points_y),g3_x(points_x,points_y),g4_x(points_x,points_y))
allocate(ni1_x(points_x,points_y),ni2_x(points_x,points_y),ni3_x(points_x,points_y),ni4_x(points_x,points_y))
allocate(s1(points_x,points_y),s2(points_x,points_y),s3(points_x,points_y),s4(points_x,points_y))
allocate(gtil1_x(points_x,points_y),gtil2_x(points_x,points_y),gtil3_x(points_x,points_y),gtil4_x(points_x,points_y))
allocate(gama_harten1_x(points_x,points_y),gama_harten2_x(points_x,points_y))
allocate(gama_harten3_x(points_x,points_y),gama_harten4_x(points_x,points_y))
!
! calcular medias de roe
!
    do j =1, points_y
        do i = 1, points_x  
            u_vetor(i,j) = Q2(i,j)/Q1(i,j)
            v_vetor(i,j) = Q3(i,j)/Q1(i,j)
            p_vetor(i,j) = (gama -1.0d0)*(Q4(i,j) - Q1(i,j)*(u_vetor(i,j)**2.0d0+v_vetor(i,j)**2.0d0)/2.0d0)
            c_vetor(i,j) = (gama*p_vetor(i,j)/Q1(i,j))**0.5d0
            H_vetor(i,j) = (Q4(i,j) + p_vetor(i,j))/Q1(i,j)
        enddo
    enddo
    !
    !
    do j = 1, points_y
        do i = 1, points_x -1
            u_roe(i,j) = (Q1(i,j)**0.5d0*u_vetor(i,j) + Q1(i+1,j)**0.5d0*u_vetor(i+1,j))/(Q1(i,j)**0.5d0 + Q1(i+1,j)**0.5d0)
            v_roe(i,j) = (Q1(i,j)**0.5d0*v_vetor(i,j) + Q1(i+1,j)**0.5d0*v_vetor(i+1,j))/(Q1(i,j)**0.5d0 + Q1(i+1,j)**0.5d0)
            H_roe(i,j) = (Q1(i,j)**0.5d0*H_vetor(i,j) + Q1(i+1,j)**0.5d0*H_vetor(i+1,j))/(Q1(i,j)**0.5d0 + Q1(i+1,j)**0.5d0)
            c_roe(i,j) = ((gama - 1.0d0)*(H_roe(i,j) - 0.50d0*(u_roe(i,j)**2.0d0 + v_roe(i,j)**2.0d0) ))**0.5d0 !WRONG (o certo é 0.5 * rho *(u² +v²) )
        enddo
    enddo
    !
    !
    do i = 1, points_x - 1
        do j = 1, points_y
            !
            !
            eigen1 = u_roe(i,j) - c_roe(i,j)
            eigen2 = u_roe(i,j)
            eigen3 = u_roe(i,j) + c_roe(i,j)
            eigen4 = u_roe(i,j)
            !
            !
            ni1_x(i,j) = (deltat/Deltax)*(u_roe(i,j) - c_roe(i,j))
            ni2_x(i,j) = (deltat/Deltax)*u_roe(i,j)
            ni3_x(i,j) = (deltat/Deltax)*(u_roe(i,j) + c_roe(i,j))
            ni4_x(i,j) = (deltat/Deltax)*u_roe(i,j)
            !
            !
            R1(1) = 1.0d0 
            R1(2) = u_roe(i,j) - c_roe(i,j)
            R1(3) = v_roe(i,j)
            R1(4) = H_roe(i,j) - u_roe(i,j)*c_roe(i,j)
            !
            !
            R2(1) = 1.0d0 
            R2(2) = u_roe(i,j)
            R2(3) = v_roe(i,j)
            R2(4) = (u_roe(i,j)**2.0d0 + v_roe(i,j)**2.0d0)/2.0d0
            !
            !
            R3(1) = 1.0d0 
            R3(2) = u_roe(i,j) + c_roe(i,j)
            R3(3) = v_roe(i,j)
            R3(4) = H_roe(i,j) + u_roe(i,j)*c_roe(i,j) 
            !
            !
            R4(1) = 0.0d0
            R4(2) = 0.0d0
            R4(3) = 1.0d0
            R4(4) = v_roe(i,j)
            !
            !
            delta_rho  = Q1(i+1,j) - Q1(i,j) 
            delta_e    = Q4(i+1,j) - Q4(i,j)
            delta_rhou = Q2(i+1,j) - Q2(i,j)
            delta_rhov = Q3(i+1,j) - Q3(i,j)
            C1 =  ((gama -1.0d0)/c_roe(i,j)**2.0d0)*(delta_e + (u_roe(i,j)**2.0d0 & 
                   + v_roe(i,j)**2.0d0)*delta_rho/2.0d0 - u_roe(i,j)*delta_rhou - v_roe(i,j)*delta_rhov )
            C2 = (delta_rhou - u_roe(i,j)*delta_rho)/c_roe(i,j)  
            !
            !
            alfa1 =  0.5d0*(C1 - C2)
            alfa2 =  delta_rho - C1 
            alfa3 =  0.5d0*(C1 + C2)
            alfa4 =  delta_rhov - v_roe(i,j)*delta_rho
            !
            !
            gtil1_x(i,j) = 0.5d0*(q_z(ni1_x(i,j),eps1) - ni1_x(i,j)**2.0d0)*alfa1
            gtil2_x(i,j) = 0.5d0*(q_z(ni2_x(i,j),eps2) - ni2_x(i,j)**2.0d0)*alfa2
            gtil3_x(i,j) = 0.5d0*(q_z(ni3_x(i,j),eps3) - ni3_x(i,j)**2.0d0)*alfa3
            gtil4_x(i,j) = 0.5d0*(q_z(ni4_x(i,j),eps4) - ni4_x(i,j)**2.0d0)*alfa4
            !
            if(i /= 1)then
                !
                !
                ni1_x(i-1,j) = (deltat/Deltax)*(u_roe(i-1,j) - c_roe(i-1,j))
                ni2_x(i-1,j) = (deltat/Deltax)*u_roe(i-1,j)
                ni3_x(i-1,j) = (deltat/Deltax)*(u_roe(i-1,j) + c_roe(i-1,j))
                ni4_x(i-1,j) = (deltat/Deltax)*u_roe(i-1,j)
                !
                gtil1_x(i-1,j) = 0.5d0*(q_z(ni1_x(i-1,j),eps1) - ni1_x(i-1,j)**2.0d0)*alfa1
                gtil2_x(i-1,j) = 0.5d0*(q_z(ni2_x(i-1,j),eps2) - ni2_x(i-1,j)**2.0d0)*alfa2
                gtil3_x(i-1,j) = 0.5d0*(q_z(ni3_x(i-1,j),eps3) - ni3_x(i-1,j)**2.0d0)*alfa3
                gtil4_x(i-1,j) = 0.5d0*(q_z(ni4_x(i-1,j),eps4) - ni4_x(i-1,j)**2.0d0)*alfa4
                !
                !
            endif
            if (abs(gtil1_x(i,j)) > 1e-6 )then
                s1(i,j) = gtil1_x(i,j)/abs(gtil1_x(i,j))
            endif
            if (abs(gtil2_x(i,j)) > 1e-6 )then
                s2(i,j) = gtil2_x(i,j)/abs(gtil2_x(i,j))
            endif
            if (abs(gtil3_x(i,j)) > 1e-6 )then
                s3(i,j) = gtil3_x(i,j)/abs(gtil3_x(i,j))
            endif
            if (abs(gtil4_x(i,j)) > 1e-6 )then
                s4(i,j) = gtil4_x(i,j)/abs(gtil4_x(i,j))
            endif
            !
            !
            if(i == 1)then
                g1_x(i,j)   = s1(i,j)  *max(0.0d0,min(abs(gtil1_x(i,j)),0.0d0))
                g2_x(i,j)   = s2(i,j)  *max(0.0d0,min(abs(gtil2_x(i,j)),0.0d0))
                g3_x(i,j)   = s3(i,j)  *max(0.0d0,min(abs(gtil3_x(i,j)),0.0d0))
                g4_x(i,j)   = s4(i,j)  *max(0.0d0,min(abs(gtil4_x(i,j)),0.0d0))
                else
                !
                g1_x(i,j)   = s1(i,j)  *max(0.0d0,min(abs(gtil1_x(i,j)),gtil1_x(i-1,j)*s1(i,j)))
                g2_x(i,j)   = s2(i,j)  *max(0.0d0,min(abs(gtil2_x(i,j)),gtil2_x(i-1,j)*s2(i,j)))
                g3_x(i,j)   = s3(i,j)  *max(0.0d0,min(abs(gtil3_x(i,j)),gtil3_x(i-1,j)*s3(i,j)))
                g4_x(i,j)   = s4(i,j)  *max(0.0d0,min(abs(gtil4_x(i,j)),gtil4_x(i-1,j)*s4(i,j)))
                !
            endif
            !
            g1_x(i+1,j)   = s1(i+1,j)  *max(0.0d0,min(abs(gtil1_x(i+1,j)),gtil1_x(i,j)*s1(i+1,j)))
            g2_x(i+1,j)   = s2(i+1,j)  *max(0.0d0,min(abs(gtil2_x(i+1,j)),gtil2_x(i,j)*s2(i+1,j)))
            g3_x(i+1,j)   = s3(i+1,j)  *max(0.0d0,min(abs(gtil3_x(i+1,j)),gtil3_x(i,j)*s3(i+1,j)))
            g4_x(i+1,j)   = s4(i+1,j)  *max(0.0d0,min(abs(gtil4_x(i+1,j)),gtil4_x(i,j)*s4(i+1,j)))
            !
            !
            if (alfa1 >= 1e-6) then
                gama_harten1_x(i,j) = (g1_x(i+1,j) -g1_x(i,j))/alfa1
            else
                gama_harten1_x(i,j) = 0.0d0
            endif
            !
            if (alfa2 >= 1e-6) then
                gama_harten2_x(i,j) = (g2_x(i+1,j) -g2_x(i,j))/alfa2
            else
                gama_harten2_x(i,j) = 0.0d0
            endif
            !
            if (alfa1 >= 1e-6) then
                gama_harten3_x(i,j) = (g3_x(i+1,j) -g3_x(i,j))/alfa3
            else
                gama_harten3_x(i,j) = 0.0d0
            endif
            !
            if (alfa4 >= 1e-6) then
                gama_harten4_x(i,j) = (g4_x(i+1,j) -g4_x(i,j))/alfa4
            else
                gama_harten4_x(i,j) = 0.0d0
            endif
            !
            !
            arg1 = ni1_x(i,j) + gama_harten1_x(i,j)
            arg2 = ni2_x(i,j) + gama_harten2_x(i,j)
            arg3 = ni3_x(i,j) + gama_harten3_x(i,j)
            arg4 = ni4_x(i,j) + gama_harten4_x(i,j)
            !
            !
            psi(1) = q_z(arg1,eps1)
            psi(2) = q_z(arg2,eps2)
            psi(3) = q_z(arg3,eps3)
            psi(4) = q_z(arg4,eps4)
            !
            ! eh o harten de primeira ordem
            !
                                                    
            E_Harten1(i,j) = 0.50d0*(E_flux1(i,j) + E_flux1(i+1,j) + (Deltax/deltat)*R1(1)*(g1_x(i,j) + g1_x(i+1,j) - psi(1)*alfa1) &
                                                                   + (Deltax/deltat)*R2(1)*(g2_x(i,j) + g2_x(i+1,j) - psi(2)*alfa2) &
                                                                   + (Deltax/deltat)*R3(1)*(g3_x(i,j) + g3_x(i+1,j) - psi(3)*alfa3) &
                                                                   + (Deltax/deltat)*R4(1)*(g4_x(i,j) + g4_x(i+1,j) - psi(4)*alfa4))
            E_Harten2(i,j) = 0.50d0*(E_flux2(i,j) + E_flux2(i+1,j) + (Deltax/deltat)*R1(2)*(g1_x(i,j) + g1_x(i+1,j) - psi(1)*alfa1) & 
                                                                   + (Deltax/deltat)*R2(2)*(g2_x(i,j) + g2_x(i+1,j) - psi(2)*alfa2) & 
                                                                   + (Deltax/deltat)*R3(2)*(g3_x(i,j) + g3_x(i+1,j) - psi(3)*alfa3) & 
                                                                   + (Deltax/deltat)*R4(2)*(g4_x(i,j) + g4_x(i+1,j) - psi(4)*alfa4))
            E_Harten3(i,j) = 0.50d0*(E_flux3(i,j) + E_flux3(i+1,j) + (Deltax/deltat)*R1(3)*(g1_x(i,j) + g1_x(i+1,j) - psi(1)*alfa1) &
                                                                   + (Deltax/deltat)*R2(3)*(g2_x(i,j) + g2_x(i+1,j) - psi(2)*alfa2) & 
                                                                   + (Deltax/deltat)*R3(3)*(g3_x(i,j) + g3_x(i+1,j) - psi(3)*alfa3) &
                                                                   + (Deltax/deltat)*R4(3)*(g4_x(i,j) + g4_x(i+1,j) - psi(4)*alfa4))
            E_Harten4(i,j) = 0.50d0*(E_flux4(i,j) + E_flux4(i+1,j) + (Deltax/deltat)*R1(4)*(g1_x(i,j) + g1_x(i+1,j) - psi(1)*alfa1) &
                                                                   + (Deltax/deltat)*R2(4)*(g2_x(i,j) + g2_x(i+1,j) - psi(2)*alfa2) &
                                                                   + (Deltax/deltat)*R3(4)*(g3_x(i,j) + g3_x(i+1,j) - psi(3)*alfa3) &
                                                                   + (Deltax/deltat)*R4(4)*(g4_x(i,j) + g4_x(i+1,j) - psi(4)*alfa4))
        enddo
    enddo
!
!
deallocate(g1_x,g2_x,g3_x,g4_x)
deallocate(ni1_x,ni2_x,ni3_x,ni4_x)
deallocate(s1,s2,s3,s4)
deallocate(gtil1_x,gtil2_x,gtil3_x,gtil4_x)
deallocate(gama_harten1_x,gama_harten2_x,gama_harten3_x,gama_harten4_x)
!
!
end subroutine fluxHarten_x
!
!**********fluxos de Harten na direcao y**********
!
subroutine fluxHarten_y
use vars
implicit none
!
!
real(8), dimension(:,:), allocatable            :: g1_y,g2_y,g3_y,g4_y
real(8), dimension(:,:), allocatable            :: gtil1_y,gtil2_y,gtil3_y,gtil4_y
real(8), dimension(:,:), allocatable            :: s1,s2,s3,s4
real(8), dimension(:,:), allocatable            :: gama_harten1_y,gama_harten2_y,gama_harten3_y,gama_harten4_y
real(8), dimension(:,:), allocatable            :: ni1_y,ni2_y,ni3_y,ni4_y 
!
! eps
!
    eps1 = 0.01d0
    eps2 = 0.0d0
    eps3 = 0.0d0
    eps4 = 0.0d0
!
!
allocate(g1_y(points_x,points_y),g2_y(points_x,points_y),g3_y(points_x,points_y),g4_y(points_x,points_y))
allocate(ni1_y(points_x,points_y),ni2_y(points_x,points_y),ni3_y(points_x,points_y),ni4_y(points_x,points_y))
allocate(s1(points_x,points_y),s2(points_x,points_y),s3(points_x,points_y),s4(points_x,points_y))
allocate(gtil1_y(points_x,points_y),gtil2_y(points_x,points_y),gtil3_y(points_x,points_y),gtil4_y(points_x,points_y))
allocate(gama_harten1_y(points_x,points_y),gama_harten2_y(points_x,points_y))
allocate(gama_harten3_y(points_x,points_y),gama_harten4_y(points_x,points_y))
!
! calcular medias de roe
!
    do i =1, points_x
        do j = 1, points_y  
            u_vetor(i,j) = Q2(i,j)/Q1(i,j)
            v_vetor(i,j)=  Q3(i,j)/Q1(i,j)
            p_vetor(i,j) = (gama -1.0d0)*(Q4(i,j) - Q1(i,j)*(u_vetor(i,j)**2.0d0+v_vetor(i,j)**2.0d0)/2.0d0)
            c_vetor(i,j) = (gama*p_vetor(i,j)/Q1(i,j))**0.5d0
            H_vetor(i,j) = (Q4(i,j) + p_vetor(i,j))/Q1(i,j)
        enddo
    enddo
    !
    !
    do i = 1, points_x
        do j = 1, points_y -1
            u_roe(i,j) = (Q1(i,j)**0.5d0*u_vetor(i,j) + Q1(i,j+1)**0.5d0*u_vetor(i,j+1))/(Q1(i,j)**0.5d0 + Q1(i,j+1)**0.5d0)
            v_roe(i,j) = (Q1(i,j)**0.5d0*v_vetor(i,j) + Q1(i,j+1)**0.5d0*v_vetor(i,j+1))/(Q1(i,j)**0.5d0 + Q1(i,j+1)**0.5d0)
            H_roe(i,j) = (Q1(i,j)**0.5d0*H_vetor(i,j) + Q1(i,j+1)**0.5d0*H_vetor(i,j+1))/(Q1(i,j)**0.5d0 + Q1(i,j+1)**0.5d0)
            c_roe(i,j) = ((gama - 1.0d0)*(H_roe(i,j) - 0.50d0*(u_roe(i,j)**2.0d0 + v_roe(i,j)**2.0d0)))**0.5d0 !WRONG (ARROOOOOW de novo)
        enddo
    enddo
    !
    !
    do i = 1, points_x  
        do j = 1, points_y - 1 
            !
            !
            eigen1 = v_roe(i,j) - c_roe(i,j)
            eigen2 = v_roe(i,j)
            eigen3 = v_roe(i,j) + c_roe(i,j)
            eigen4 = v_roe(i,j)         
            !
            !
            ni1_y(i,j) = (deltat/Deltay)*(v_roe(i,j) - c_roe(i,j))
            ni2_y(i,j) = (deltat/Deltay)*v_roe(i,j)
            ni3_y(i,j) = (deltat/Deltay)*(v_roe(i,j) + c_roe(i,j))
            ni4_y(i,j) = (deltat/Deltay)*v_roe(i,j)
            !
            !
            R1(1) = 1.0d0 
            R1(2) = u_roe(i,j) 
            R1(3) = v_roe(i,j) - c_roe(i,j)
            R1(4) = H_roe(i,j) - v_roe(i,j)*c_roe(i,j)
            !
            !
            R2(1) = 1.0d0 
            R2(2) = u_roe(i,j)
            R2(3) = v_roe(i,j)
            R2(4) = (u_roe(i,j)**2.0d0 + v_roe(i,j)**2.0d0)/2.0d0
            !
            !
            R3(1) = 1.0d0
            R3(2) = u_roe(i,j) 
            R3(3) = v_roe(i,j) + c_roe(i,j)
            R3(4) = H_roe(i,j) + v_roe(i,j)*c_roe(i,j)
            !
            !
            R4(1) = 0.0d0 
            R4(2) = 1.0d0
            R4(3) = 0.0d0
            R4(4) = u_roe(i,j)
            !
            !
            delta_rho  = Q1(i,j+1) - Q1(i,j) 
            delta_e    = Q4(i,j+1) - Q4(i,j)
            delta_rhou = Q2(i,j+1) - Q2(i,j)
            delta_rhov = Q3(i,j+1) - Q3(i,j)
            C1 =  ((gama -1.0d0)/c_roe(i,j)**2.0d0)*(delta_e + (u_roe(i,j)**2.0d0 & 
                  + v_roe(i,j)**2.0d0)*delta_rho/2.0d0 - u_roe(i,j)*delta_rhou - v_roe(i,j)*delta_rhov )
            C2 = (delta_rhov - v_roe(i,j)*delta_rho)/c_roe(i,j)  
            !
            !
            alfa1 =  0.5d0*(C1 - C2)
            alfa2 =  delta_rho - C1 
            alfa3 =  0.5d0*(C1 + C2)
            alfa4 =  delta_rhou - u_roe(i,j)*delta_rho
            !
            !
            gtil1_y(i,j) = 0.5d0*(q_z(ni1_y(i,j),eps1) - ni1_y(i,j)**2.0d0)*alfa1
            gtil2_y(i,j) = 0.5d0*(q_z(ni2_y(i,j),eps2) - ni2_y(i,j)**2.0d0)*alfa2
            gtil3_y(i,j) = 0.5d0*(q_z(ni3_y(i,j),eps3) - ni3_y(i,j)**2.0d0)*alfa3
            gtil4_y(i,j) = 0.5d0*(q_z(ni4_y(i,j),eps4) - ni4_y(i,j)**2.0d0)*alfa4
            !
            if (j /= 1)then
                ni1_y(i,j-1) = (deltat/Deltay)*(v_roe(i,j-1) - c_roe(i,j-1))
                ni2_y(i,j-1) = (deltat/Deltay)*v_roe(i,j-1)
                ni3_y(i,j-1) = (deltat/Deltay)*(v_roe(i,j-1) + c_roe(i,j-1))
                ni4_y(i,j-1) = (deltat/Deltay)*v_roe(i,j-1)
                !
                gtil1_y(i,j-1) = 0.5d0*(q_z(ni1_y(i,j-1),eps1) - ni1_y(i,j)**2.0d0)*alfa1
                gtil2_y(i,j-1) = 0.5d0*(q_z(ni2_y(i,j-1),eps2) - ni2_y(i,j)**2.0d0)*alfa2
                gtil3_y(i,j-1) = 0.5d0*(q_z(ni3_y(i,j-1),eps3) - ni3_y(i,j)**2.0d0)*alfa3
                gtil4_y(i,j-1) = 0.5d0*(q_z(ni4_y(i,j-1),eps4) - ni4_y(i,j)**2.0d0)*alfa4
            endif
            !
            !
            if (abs(gtil1_y(i,j)) > 1e-6 )then
                s1(i,j) = gtil1_y(i,j)/abs(gtil1_y(i,j))
            endif
            if (abs(gtil2_y(i,j)) > 1e-6 )then
                s2(i,j) = gtil2_y(i,j)/abs(gtil2_y(i,j))
            endif
            if (abs(gtil3_y(i,j)) > 1e-6 )then
                s3(i,j) = gtil3_y(i,j)/abs(gtil3_y(i,j))
            endif
            if (abs(gtil4_y(i,j)) > 1e-6 )then
                s4(i,j) = gtil4_y(i,j)/abs(gtil4_y(i,j))
            endif
            !
            !
            if (j == 1)then
                g1_y(i,j)   = s1(i,j)  *max(0.0d0,min(abs(gtil1_y(i,j)),0.0d0))
                g2_y(i,j)   = s2(i,j)  *max(0.0d0,min(abs(gtil2_y(i,j)),0.0d0))
                g3_y(i,j)   = s3(i,j)  *max(0.0d0,min(abs(gtil3_y(i,j)),0.0d0))
                g4_y(i,j)   = s4(i,j)  *max(0.0d0,min(abs(gtil4_y(i,j)),0.0d0))
                else
                g1_y(i,j)   = s1(i,j)  *max(0.0d0,min(abs(gtil1_y(i,j)),gtil1_y(i,j-1)*s1(i,j)))
                g2_y(i,j)   = s2(i,j)  *max(0.0d0,min(abs(gtil2_y(i,j)),gtil2_y(i,j-1)*s2(i,j)))
                g3_y(i,j)   = s3(i,j)  *max(0.0d0,min(abs(gtil3_y(i,j)),gtil3_y(i,j-1)*s3(i,j)))
                g4_y(i,j)   = s4(i,j)  *max(0.0d0,min(abs(gtil4_y(i,j)),gtil4_y(i,j-1)*s4(i,j)))    
            endif
            !
            g1_y(i,j+1)   = s1(i,j+1)  *max(0.0d0,min(abs(gtil1_y(i,j+1)),gtil1_y(i,j)*s1(i,j+1)))
            g2_y(i,j+1)   = s2(i,j+1)  *max(0.0d0,min(abs(gtil2_y(i,j+1)),gtil2_y(i,j)*s2(i,j+1)))
            g3_y(i,j+1)   = s3(i,j+1)  *max(0.0d0,min(abs(gtil3_y(i,j+1)),gtil3_y(i,j)*s3(i,j+1)))
            g4_y(i,j+1)   = s4(i,j+1)  *max(0.0d0,min(abs(gtil4_y(i,j+1)),gtil4_y(i,j)*s4(i,j+1)))
            !
            !
            if (alfa1 >= 1e-6) then
                gama_harten1_y(i,j) = (g1_y(i,j+1) -g1_y(i,j))/alfa1
            else
                gama_harten1_y(i,j) = 0.0d0
            endif
            !
            if (alfa2 >= 1e-6) then
                gama_harten2_y(i,j) = (g2_y(i,j+1) -g2_y(i,j))/alfa2
            else
                gama_harten2_y(i,j) = 0.0d0
            endif
            !
            if (alfa1 >= 1e-6) then
                gama_harten3_y(i,j) = (g3_y(i,j+1) -g3_y(i,j))/alfa3
            else
                gama_harten3_y(i,j) = 0.0d0
            endif
            !
            if (alfa4 >= 1e-6) then
                gama_harten4_y(i,j) = (g4_y(i,j+1) -g4_y(i,j))/alfa4
            else
                gama_harten4_y(i,j) = 0.0d0
            endif
            !
            !
            arg1 = ni1_y(i,j) + gama_harten1_y(i,j)
            arg2 = ni2_y(i,j) + gama_harten2_y(i,j)
            arg3 = ni3_y(i,j) + gama_harten3_y(i,j)
            arg4 = ni4_y(i,j) + gama_harten4_y(i,j)
            !
            !
            psi(1) = q_z(arg1,eps1)
            psi(2) = q_z(arg2,eps2)
            psi(3) = q_z(arg3,eps3)
            psi(4) = q_z(arg4,eps4)
            !
            ! 
            !
                                                    
            F_Harten1(i,j) = 0.50d0*(F_flux1(i,j) + F_flux1(i,j+1) + (Deltay/deltat)*R1(1)*(g1_y(i,j) + g1_y(i,j+1) - psi(1)*alfa1) & 
                                                                   + (Deltay/deltat)*R2(1)*(g2_y(i,j) + g2_y(i,j+1) - psi(2)*alfa2) &
                                                                   + (Deltay/deltat)*R3(1)*(g3_y(i,j) + g3_y(i,j+1) - psi(3)*alfa3) &
                                                                   + (Deltay/deltat)*R4(1)*(g4_y(i,j) + g4_y(i,j+1) - psi(4)*alfa4))
            F_Harten2(i,j) = 0.50d0*(F_flux2(i,j) + F_flux2(i,j+1) + (Deltay/deltat)*R1(2)*(g1_y(i,j) + g1_y(i,j+1) - psi(1)*alfa1) &
                                                                   + (Deltay/deltat)*R2(2)*(g2_y(i,j) + g2_y(i,j+1) - psi(2)*alfa2) &
                                                                   + (Deltay/deltat)*R3(2)*(g3_y(i,j) + g3_y(i,j+1) - psi(3)*alfa3) &
                                                                   + (Deltay/deltat)*R4(2)*(g4_y(i,j) + g4_y(i,j+1) - psi(4)*alfa4))
            F_Harten3(i,j) = 0.50d0*(F_flux3(i,j) + F_flux3(i,j+1) + (Deltay/deltat)*R1(3)*(g1_y(i,j) + g1_y(i,j+1) - psi(1)*alfa1) &
                                                                   + (Deltay/deltat)*R2(3)*(g2_y(i,j) + g2_y(i,j+1) - psi(2)*alfa2) &
                                                                   + (Deltay/deltat)*R3(3)*(g3_y(i,j) + g3_y(i,j+1) - psi(3)*alfa3) &
                                                                   + (Deltay/deltat)*R4(3)*(g4_y(i,j) + g4_y(i,j+1) - psi(4)*alfa4))
            F_Harten4(i,j) = 0.50d0*(F_flux4(i,j) + F_flux4(i,j+1) + (Deltay/deltat)*R1(4)*(g1_y(i,j) + g1_y(i,j+1) - psi(1)*alfa1) &
                                                                   + (Deltay/deltat)*R2(4)*(g2_y(i,j) + g2_y(i,j+1) - psi(2)*alfa2) &
                                                                   + (Deltay/deltat)*R3(4)*(g3_y(i,j) + g3_y(i,j+1) - psi(3)*alfa3) &
                                                                   + (Deltay/deltat)*R4(4)*(g4_y(i,j) + g4_y(i,j+1) - psi(4)*alfa4))
        enddo
    enddo
!
!
deallocate(g1_y,g2_y,g3_y,g4_y)
deallocate(ni1_y,ni2_y,ni3_y,ni4_y)
deallocate(s1,s2,s3,s4)
deallocate(gtil1_y,gtil2_y,gtil3_y,gtil4_y)
deallocate(gama_harten1_y,gama_harten2_y,gama_harten3_y,gama_harten4_y)
!
!
end subroutine fluxHarten_y
!
!
!
!subroutine Harten_implicit
!use vars
!implicit none

!end subroutine Harten_implicit
!
!
!
subroutine convergence
use vars
implicit none
!
!
!
    do i = 2, points_x - 1
        do j = 2, points_y -1
            !
            !
            residue1(i,j) = (1.0d0/Deltax)*(E_Harten1(i,j) - E_Harten1(i-1,j)) + (1.0d0/Deltay)*(F_Harten1(i,j) - F_Harten1(i,j-1)) 
            residue2(i,j) = (1.0d0/Deltax)*(E_Harten2(i,j) - E_Harten2(i-1,j)) + (1.0d0/Deltay)*(F_Harten2(i,j) - F_Harten2(i,j-1)) 
            residue3(i,j) = (1.0d0/Deltax)*(E_Harten3(i,j) - E_Harten3(i-1,j)) + (1.0d0/Deltay)*(F_Harten3(i,j) - F_Harten3(i,j-1)) 
            residue4(i,j) = (1.0d0/Deltax)*(E_Harten4(i,j) - E_Harten4(i-1,j)) + (1.0d0/Deltay)*(F_Harten4(i,j) - F_Harten4(i,j-1)) 
            !
            !
            residue = max(residue1(i,j),residue2(i,j),residue3(i,j),residue4(i,j))
            !
        enddo
    enddo
!
!
!
end subroutine convergence
!
!**********HARTEN**********
!
subroutine Harten
use vars
implicit none
!
!
real(8),dimension(:,:), allocatable         :: u_max, v_max
real(8),dimension(:,:), allocatable         :: p_max, a_max,CFL_vetor
!
!
!
allocate(u_max(points_x,points_y),v_max(points_x,points_y))
allocate(p_max(points_x,points_y),a_max(points_x,points_y),CFL_vetor(points_x,points_y))
!
!
    do i = 2, points_x -1
        do j = 2, points_y - 1
            Q1(i,j) = Q1(i,j) - deltat*residue1(i,j) 
            Q2(i,j) = Q2(i,j) - deltat*residue2(i,j)
            Q3(i,j) = Q3(i,j) - deltat*residue3(i,j)
            Q4(i,j) = Q4(i,j) - deltat*residue4(i,j)
        enddo
    enddo
!
! calculo de CFL
!    
	do i = 1,points_x
		do j = 1,points_y
			u_max(i,j) = (Q2(i,j)/Q1(i,j))
			v_max(i,j) = (Q3(i,j)/Q1(i,j))
			p_max(i,j) = (gama -1.0d0)*(Q4(i,j) - Q1(i,j)*(u_max(i,j)**2.0d0+v_max(i,j)**2.0d0)/2.0d0)
			a_max(i,j) = sqrt(gama*p_max(i,j)/Q1(i,j))
			CFL_vetor(i,j) = (deltat/Deltax)*(u_max(i,j)+a_max(i,j))
		enddo
	enddo
!
!
	CFL = maxval(CFL_vetor)
!
!
deallocate(u_max,v_max)
deallocate(p_max,a_max,CFL_vetor)
!
!
end subroutine Harten
!
!**********Arquivos de saida**********
!
subroutine output
use vars
implicit none
!
!
real(8),dimension(:,:), allocatable         :: u_out, v_out
real(8),dimension(:,:), allocatable         :: p_out
!
!
!       open(10,file='proj4')
!
!       do i = 1, malha + 2
!           WRITE(10,*) x(i),Q(i,2)/Q(i,1)
!       end do
!
allocate(u_out(points_x,points_y),v_out(points_x,points_y))
allocate(p_out(points_x,points_y))
!
!
do j = 1, points_y
    do i = 1,points_x
        u_out(i,j) = Q2(i,j)/Q1(i,j)
        v_out(i,j) = Q3(i,j)/Q1(i,j)
    enddo
enddo
!
!
do j = 1, points_y
    do i = 1,points_x
        p_out(i,j) = (gama -1.0d0)*(Q4(i,j) - Q1(i,j)*(u_out(i,j)**2.0d0+v_out(i,j)**2.0d0)/2.0d0)
    enddo
enddo
!
!gravar arquivo para o Tecplot.
!
!open(4,file='malha.dat')
!!
!write(4,*)' TITLE="malha do projeto" '
!write(4,*)' VARIABLES = "x", "y" ,"u", "v", "rho" "p" "e" "T/T0"'
!write(4,*)' ZONE I=',points_x,' J=',points_y,',DATAPACKING=POINT'
!! 
!do j = 1,points_y
!   do i = 1,points_x
!      WRITE(4,'(10es11.3e2)') meshx(i,j),meshy(i,j),u_out(i,j),v_out(i,j),Q1(i,j),p_out(i,j),Q4(i,j)
!   enddo
!enddo
!!
!!
!close(4)
!
open(10,file='resultados_numericos')
        j = 40
        do i = 1, points_x
                WRITE(10,*) mach_inlet,p_percent,shock_angle,gama
        enddo
!
close(10)
!
!
deallocate(u_out,v_out)
deallocate(p_out)
!
!
end subroutine output
!
!**********Arquivos de saida**********
!
subroutine transient
use vars
implicit none
!
!
real(8),dimension(:,:), allocatable         :: u_out_trans, v_out_trans, T
real(8),dimension(:,:), allocatable         :: p_out_trans, temp_estag, M_final
!
!
!       open(10,file='proj4')
!
!       do i = 1, malha + 2
!           WRITE(10,*) x(i),Q(i,2)/Q(i,1)
!       end do
!
allocate(u_out_trans(points_x,points_y),v_out_trans(points_x,points_y),T(points_x,points_y))
allocate(p_out_trans(points_x,points_y),temp_estag(points_x,points_y),M_final(points_x,points_y))
!
!
!
do i = 1,points_x
    do j = 1, points_y
        u_out_trans(i,j) = Q2(i,j)/Q1(i,j)
        v_out_trans(i,j) = Q3(i,j)/Q1(i,j)
    enddo
enddo
!
!
do i = 1,points_x
    do j = 1, points_y
        M_final(i,j) = sqrt(u_out_trans(i,j)**2.0d0 + v_out_trans(i,j)**2.0d0)
    enddo
enddo
!
!
do i = 1, points_x
    do j = 1, points_y
        temp_estag(i,j) = 1.0d0 + ((gama-1.0d0)/2.0d0)*(M_final(i,j)**2.0d0)
        T(i,j) = 288.d0*temp_estag(i,j)
    enddo
enddo
!
!
do i = 1,points_x
    do j = 1, points_y
        p_out_trans(i,j) = (gama -1.0d0)*(Q4(i,j) - Q1(i,j)*(u_out_trans(i,j)**2.0d0+v_out_trans(i,j)**2.0d0)/2.0d0)
    enddo
enddo
!
!gravar arquivo para o Tecplot.
!
if(tec_sol == 1) then
    !
    open(4,file='1.dat')
    !
    write(4,*)' TITLE="malha do projeto" '
    write(4,*)' VARIABLES = "x", "y" ,"u", "v", "rho" "p" "e" "T/T0" "T" "Mach" '
    write(4,*)' ZONE I=',points_x,' J=',points_y,',DATAPACKING=POINT'
    !
    do j = 1,points_y
	   do i = 1,points_x
          WRITE(4,'(10es11.3e2)') meshx(i,j),meshy(i,j),u_out_trans(i,j),v_out_trans(i,j),Q1(i,j),p_out_trans(i,j),Q4(i,j),temp_estag(i,j), T(i,j),M_final(i,j)
       enddo
    enddo
    !
    !
close(4)
else if(tec_sol == 2) then
    !
    open(4,file='2.dat')
    !
    write(4,*)' TITLE="malha do projeto" '
    write(4,*)' VARIABLES = "x", "y" ,"u", "v", "rho" "p" "e" "T/T0" "T" "Mach" '
    write(4,*)' ZONE I=',points_x,' J=',points_y,',DATAPACKING=POINT'
    !
    do j = 1,points_y
	   do i = 1,points_x
          WRITE(4,'(10es11.3e2)') meshx(i,j),meshy(i,j),u_out_trans(i,j),v_out_trans(i,j),Q1(i,j),p_out_trans(i,j),Q4(i,j),temp_estag(i,j), T(i,j),M_final(i,j)
       enddo
    enddo
    !
    !
close(4)
else if(tec_sol == 3) then
    !
    open(4,file='3.dat')
    !
    write(4,*)' TITLE="malha do projeto" '
    write(4,*)' VARIABLES = "x", "y" ,"u", "v", "rho" "p" "e" "T/T0" "T" "Mach" '
    write(4,*)' ZONE I=',points_x,' J=',points_y,',DATAPACKING=POINT'
    !
    do j = 1,points_y
	   do i = 1,points_x
          WRITE(4,'(10es11.3e2)') meshx(i,j),meshy(i,j),u_out_trans(i,j),v_out_trans(i,j),Q1(i,j),p_out_trans(i,j),Q4(i,j),temp_estag(i,j), T(i,j),M_final(i,j)
       enddo
    enddo
    !
    !
close(4)
else if(tec_sol == 4) then
    !
    open(4,file='4.dat')
    !
    write(4,*)' TITLE="malha do projeto" '
    write(4,*)' VARIABLES = "x", "y" ,"u", "v", "rho" "p" "e" "T/T0" "T" "Mach" '
    write(4,*)' ZONE I=',points_x,' J=',points_y,',DATAPACKING=POINT'
    !
    do j = 1,points_y
	   do i = 1,points_x
          WRITE(4,'(10es11.3e2)') meshx(i,j),meshy(i,j),u_out_trans(i,j),v_out_trans(i,j),Q1(i,j),p_out_trans(i,j),Q4(i,j),temp_estag(i,j), T(i,j),M_final(i,j)
       enddo
    enddo
    !
    !
close(4)
else if(tec_sol == 5) then
    open(4,file='5.dat')
    !
    write(4,*)' TITLE="malha do projeto" '
    write(4,*)' VARIABLES = "x", "y" ,"u", "v", "rho" "p" "e" "T/T0" "T" "Mach" '
    write(4,*)' ZONE I=',points_x,' J=',points_y,',DATAPACKING=POINT'
    !
    do j = 1,points_y
	   do i = 1,points_x
          WRITE(4,'(10es11.3e2)') meshx(i,j),meshy(i,j),u_out_trans(i,j),v_out_trans(i,j),Q1(i,j),p_out_trans(i,j),Q4(i,j),temp_estag(i,j), T(i,j),M_final(i,j)
       enddo
    enddo
    !
    !
close(4)
else if(tec_sol == 6) then
    open(4,file='6.dat')
    !
    write(4,*)' TITLE="malha do projeto" '
    write(4,*)' VARIABLES = "x", "y" ,"u", "v", "rho" "p" "e" "T/T0" "T" "Mach" '
    write(4,*)' ZONE I=',points_x,' J=',points_y,',DATAPACKING=POINT'
    !
    do j = 1,points_y
	   do i = 1,points_x
          WRITE(4,'(10es11.3e2)') meshx(i,j),meshy(i,j),u_out_trans(i,j),v_out_trans(i,j),Q1(i,j),p_out_trans(i,j),Q4(i,j),temp_estag(i,j), T(i,j),M_final(i,j)
       enddo
    enddo
    !
    !
close(4)
else if(tec_sol == 7) then
    open(4,file='7.dat')
    !
    write(4,*)' TITLE="malha do projeto" '
    write(4,*)' VARIABLES = "x", "y" ,"u", "v", "rho" "p" "e" "T/T0" "T" "Mach" '
    write(4,*)' ZONE I=',points_x,' J=',points_y,',DATAPACKING=POINT'
    !
    do j = 1,points_y
	   do i = 1,points_x
          WRITE(4,'(10es11.3e2)') meshx(i,j),meshy(i,j),u_out_trans(i,j),v_out_trans(i,j),Q1(i,j),p_out_trans(i,j),Q4(i,j),temp_estag(i,j), T(i,j),M_final(i,j)
       enddo
    enddo
    !
    !
close(4)
else if(tec_sol == 8) then
    open(4,file='8.dat')
    !
    write(4,*)' TITLE="malha do projeto" '
    write(4,*)' VARIABLES = "x", "y" ,"u", "v", "rho" "p" "e" "T/T0" "T" "Mach" '
    write(4,*)' ZONE I=',points_x,' J=',points_y,',DATAPACKING=POINT'
    !
    do j = 1,points_y
	   do i = 1,points_x
          WRITE(4,'(10es11.3e2)') meshx(i,j),meshy(i,j),u_out_trans(i,j),v_out_trans(i,j),Q1(i,j),p_out_trans(i,j),Q4(i,j),temp_estag(i,j), T(i,j),M_final(i,j)
       enddo
    enddo
    !
    !
close(4)
!
!
else if(tec_sol == 9) then
    !
    open(4,file='9.dat')
    !
    write(4,*)' TITLE="malha do projeto" '
    write(4,*)' VARIABLES = "x", "y" ,"u", "v", "rho" "p" "e" "T/T0" "T" "Mach" '
    write(4,*)' ZONE I=',points_x,' J=',points_y,',DATAPACKING=POINT'
    !
    do j = 1,points_y
	   do i = 1,points_x
          WRITE(4,'(10es11.3e2)') meshx(i,j),meshy(i,j),u_out_trans(i,j),v_out_trans(i,j),Q1(i,j),p_out_trans(i,j),Q4(i,j),temp_estag(i,j), T(i,j),M_final(i,j)
       enddo
    enddo
    !
    !
close(4)
!
!
else if(tec_sol == 10) then
    !
    open(4,file='10.dat')
    !
    write(4,*)' TITLE="malha do projeto" '
    write(4,*)' VARIABLES = "x", "y" ,"u", "v", "rho" "p" "e" "T/T0" "T" "Mach" '
    write(4,*)' ZONE I=',points_x,' J=',points_y,',DATAPACKING=POINT'
    !
    do j = 1,points_y
	   do i = 1,points_x
          WRITE(4,'(10es11.3e2)') meshx(i,j),meshy(i,j),u_out_trans(i,j),v_out_trans(i,j),Q1(i,j),p_out_trans(i,j),Q4(i,j),temp_estag(i,j), T(i,j),M_final(i,j)
       enddo
    enddo
    !
    !
close(4)
!
!
else if(tec_sol == 11) then
    !
    open(4,file='11.dat')
    !
    write(4,*)' TITLE="malha do projeto" '
    write(4,*)' VARIABLES = "x", "y" ,"u", "v", "rho" "p" "e" "T/T0" "T" "Mach" '
    write(4,*)' ZONE I=',points_x,' J=',points_y,',DATAPACKING=POINT'
    !
    do j = 1,points_y
	   do i = 1,points_x
          WRITE(4,'(10es11.3e2)') meshx(i,j),meshy(i,j),u_out_trans(i,j),v_out_trans(i,j),Q1(i,j),p_out_trans(i,j),Q4(i,j),temp_estag(i,j), T(i,j),M_final(i,j)
       enddo
    enddo
    !
    !
close(4)
!
!
!
!
else if(tec_sol == 12) then
    !
    open(4,file='12.dat')
    !
    write(4,*)' TITLE="malha do projeto" '
    write(4,*)' VARIABLES = "x", "y" ,"u", "v", "rho" "p" "e" "T/T0" "T" "Mach" '
    write(4,*)' ZONE I=',points_x,' J=',points_y,',DATAPACKING=POINT'
    !
    do j = 1,points_y
	   do i = 1,points_x
          WRITE(4,'(10es11.3e2)') meshx(i,j),meshy(i,j),u_out_trans(i,j),v_out_trans(i,j),Q1(i,j),p_out_trans(i,j),Q4(i,j),temp_estag(i,j), T(i,j),M_final(i,j)
       enddo
    enddo
    !
    !
close(4)
!
!
else if(tec_sol == 13) then
    !
    open(4,file='13.dat')
    !
    write(4,*)' TITLE="malha do projeto" '
    write(4,*)' VARIABLES = "x", "y" ,"u", "v", "rho" "p" "e" "T/T0" "T" "Mach" '
    write(4,*)' ZONE I=',points_x,' J=',points_y,',DATAPACKING=POINT'
    !
    do j = 1,points_y
	   do i = 1,points_x
          WRITE(4,'(10es11.3e2)') meshx(i,j),meshy(i,j),u_out_trans(i,j),v_out_trans(i,j),Q1(i,j),p_out_trans(i,j),Q4(i,j),temp_estag(i,j), T(i,j),M_final(i,j)
       enddo
    enddo
    !
    !
close(4)
!
!
else if(tec_sol == 14) then
    !
    open(4,file='14.dat')
    !
    write(4,*)' TITLE="malha do projeto" '
    write(4,*)' VARIABLES = "x", "y" ,"u", "v", "rho" "p" "e" "T/T0" "T" "Mach" '
    write(4,*)' ZONE I=',points_x,' J=',points_y,',DATAPACKING=POINT'
    !
    do j = 1,points_y
	   do i = 1,points_x
          WRITE(4,'(10es11.3e2)') meshx(i,j),meshy(i,j),u_out_trans(i,j),v_out_trans(i,j),Q1(i,j),p_out_trans(i,j),Q4(i,j),temp_estag(i,j), T(i,j),M_final(i,j)
       enddo
    enddo
    !
    !
close(4)
!
!
else if(tec_sol == 15) then
    !
    open(4,file='15.dat')
    !
    write(4,*)' TITLE="malha do projeto" '
    write(4,*)' VARIABLES = "x", "y" ,"u", "v", "rho" "p" "e" "T/T0" "T" "Mach" '
    write(4,*)' ZONE I=',points_x,' J=',points_y,',DATAPACKING=POINT'
    !
    do j = 1,points_y
	   do i = 1,points_x
          WRITE(4,'(10es11.3e2)') meshx(i,j),meshy(i,j),u_out_trans(i,j),v_out_trans(i,j),Q1(i,j),p_out_trans(i,j),Q4(i,j),temp_estag(i,j), T(i,j),M_final(i,j)
       enddo
    enddo
    !
    !
close(4)
!
!
else if(tec_sol == 16) then
    !
    open(4,file='16.dat')
    !
    write(4,*)' TITLE="malha do projeto" '
    write(4,*)' VARIABLES = "x", "y" ,"u", "v", "rho" "p" "e" "T/T0" "T" "Mach" '
    write(4,*)' ZONE I=',points_x,' J=',points_y,',DATAPACKING=POINT'
    !
    do j = 1,points_y
	   do i = 1,points_x
          WRITE(4,'(10es11.3e2)') meshx(i,j),meshy(i,j),u_out_trans(i,j),v_out_trans(i,j),Q1(i,j),p_out_trans(i,j),Q4(i,j),temp_estag(i,j), T(i,j),M_final(i,j)
       enddo
    enddo
    !
    !
close(4)
!
!
else if(tec_sol == 17) then
    !
    open(4,file='17.dat')
    !
    write(4,*)' TITLE="malha do projeto" '
    write(4,*)' VARIABLES = "x", "y" ,"u", "v", "rho" "p" "e" "T/T0" "T" "Mach" '
    write(4,*)' ZONE I=',points_x,' J=',points_y,',DATAPACKING=POINT'
    !
    do j = 1,points_y
	   do i = 1,points_x
          WRITE(4,'(10es11.3e2)') meshx(i,j),meshy(i,j),u_out_trans(i,j),v_out_trans(i,j),Q1(i,j),p_out_trans(i,j),Q4(i,j),temp_estag(i,j), T(i,j),M_final(i,j)
       enddo
    enddo
    !
    !
close(4)
!
!
else if(tec_sol == 18) then
    !
    open(4,file='18.dat')
    !
    write(4,*)' TITLE="malha do projeto" '
    write(4,*)' VARIABLES = "x", "y" ,"u", "v", "rho" "p" "e" "T/T0" "T" "Mach" '
    write(4,*)' ZONE I=',points_x,' J=',points_y,',DATAPACKING=POINT'
    !
    do j = 1,points_y
	   do i = 1,points_x
          WRITE(4,'(10es11.3e2)') meshx(i,j),meshy(i,j),u_out_trans(i,j),v_out_trans(i,j),Q1(i,j),p_out_trans(i,j),Q4(i,j),temp_estag(i,j), T(i,j),M_final(i,j)
       enddo
    enddo
    !
    !
close(4)
!
!
else if(tec_sol == 19) then
    !
    open(4,file='19.dat')
    !
    write(4,*)' TITLE="malha do projeto" '
    write(4,*)' VARIABLES = "x", "y" ,"u", "v", "rho" "p" "e" "T/T0" "T" "Mach" '
    write(4,*)' ZONE I=',points_x,' J=',points_y,',DATAPACKING=POINT'
    !
    do j = 1,points_y
	   do i = 1,points_x
          WRITE(4,'(10es11.3e2)') meshx(i,j),meshy(i,j),u_out_trans(i,j),v_out_trans(i,j),Q1(i,j),p_out_trans(i,j),Q4(i,j),temp_estag(i,j), T(i,j),M_final(i,j)
       enddo
    enddo
    !
    !
close(4)
!
!
else 
    !
    open(4,file='20.dat')
    !
    write(4,*)' TITLE="malha do projeto" '
    write(4,*)' VARIABLES = "x", "y" ,"u", "v", "rho" "p" "e" "T/T0" "T" "Mach" '
    write(4,*)' ZONE I=',points_x,' J=',points_y,',DATAPACKING=POINT'
    !
    do j = 1,points_y
	   do i = 1,points_x
          WRITE(4,'(10es11.3e2)') meshx(i,j),meshy(i,j),u_out_trans(i,j),v_out_trans(i,j),Q1(i,j),p_out_trans(i,j),Q4(i,j),temp_estag(i,j), T(i,j),M_final(i,j)
       enddo
    enddo
    !
    !
close(4)
endif
!
!
deallocate(u_out_trans,v_out_trans,T)
deallocate(p_out_trans,M_final,temp_estag)
!
!
end subroutine transient
