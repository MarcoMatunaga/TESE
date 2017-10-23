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
integer(4)                                       :: iter, dimen, max_iter, iter_save, n_save, n_sol
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
        
        character (len=2) function convert(a)
        !
        implicit none
        !
        integer(4) :: a
        !
        if (a < 10) then
        convert = char(48+a)
        else if (a >= 10 .and. a < 20) then 
        convert = char(48+1)//char(48+mod(a,10))
        !
        !
        else if (a >= 20 .and. a < 30) then
        convert = char(48+2)//char(48+mod(a,10))
        !
        !
        else if (a >= 30 .and. a < 40) then
        convert = char(48+3)//char(48+mod(a,10))
        !
        !
        else if (a >= 40 .and. a < 50) then
        convert = char(48+4)//char(48+mod(a,10))
        !
        !
        else if (a >= 50 .and. a < 60) then
        convert = char(48+5)//char(48+mod(a,10))
        !
        !
        else if (a >= 60 .and. a < 70) then
        convert = char(48+6)//char(48+mod(a,10))
        !
        !
        else if (a >= 70 .and. a < 80) then
        convert = char(48+7)//char(48+mod(a,10))
        !
        !
        else if (a >= 80 .and. a < 90) then
        convert = char(48+8)//char(48+mod(a,10))
        !
        !
        else if (a >= 90 .and. a < 100) then
        convert = char(48+9)//char(48+mod(a,10))
        !
        !
        else
        write(*,*) 'Number os solutions saved more than 99'
        end if
        !
        !
        end function convert
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
namelist /PAR_Flow/ mach_inlet,shock_angle,gama
namelist /PAR_Time/ deltat, max_residue, max_iter, iter_save
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
! Criterio de convergencia
!
iter = 0
residue = 10
n_save = 1
!
! Delta x e y, sao respectivamente os espacamentos de malha
! em x e y 
!
Deltax = Length/DBLE(points_x)
Deltay = Height/DBLE(points_y)
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
!
do while(abs(residue) > max_residue .and. iter < max_iter)
    call fluxes
    call fluxHarten_x
    call fluxHarten_y
    call convergence
    call Harten
    call boundary_conditions
    iter = iter + 1
    if(iter == iter_save*n_save) then
         call output
         n_save = n_save + 1
    end if
    open(1,file='residue')
    write(1,*) iter, abs(residue)
enddo
    close(1)
call output_final
!
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
        do j =1, points_y
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
    do j = 1, points_y

        i = 1

        u   = mach_inlet
        v   = 0.0d0
        p   = 1.0d0/gama
        rho = 1.0d0

        Q1(i,j) = rho
        Q2(i,j) = rho*u
        Q3(i,j) = rho*v
        Q4(i,j) = p/(gama-1.0d0) + 0.5d0*rho*(u**2.0d0 + v**2.0d0)

    end do

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

    do j = 1, points_y

        i = points_x

        Q1(i,j) = Q1(i-1,j)
        Q2(i,j) = Q2(i-1,j)
        Q3(i,j) = Q3(i-1,j)
        Q4(i,j) = Q4(i-1,j)

    end do
!
! Post shock conditions 
!
        rho =  (gama + 1.0d0)*mach_inlet_n**2.0d0/(2.0d0 + (gama - 1.0d0)*mach_inlet_n**2.0d0)
        p   =  (1.0d0 + (2.0d0*gama/( gama + 1.0d0 ))*(mach_inlet_n**2.0d0 - 1.0d0) ) / gama
        a   =  (gama*p/rho)**0.5d0
        u   =  mach_2 * a * cos(turn_angle)
        v   = -mach_2 * a * sin(turn_angle)
        
    do i = 1, points_x
        
        j = points_y

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
real(8),dimension(:,:), allocatable         :: u_out_trans, v_out_trans
real(8),dimension(:,:), allocatable         :: p_out_trans, M_final, P_total, P_estag
character (len=100) :: fname = '.dat'
character (len=100) :: FileTag
!
!
!
allocate(u_out_trans(points_x,points_y),v_out_trans(points_x,points_y))
allocate(p_out_trans(points_x,points_y),M_final(points_x,points_y))
allocate(P_total(points_x,points_y),P_estag(points_x,points_y))
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
        p_out_trans(i,j) = (gama -1.0d0)*(Q4(i,j) - Q1(i,j)*(u_out_trans(i,j)**2.0d0+v_out_trans(i,j)**2.0d0)/2.0d0)
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
do i = 1,points_x
    do j = 1,points_y
       P_estag(i,j) = p_out_trans(i,j)*(1.0d0 + ((gama -1.0d0)/2.0d0)*M_final(i,j)**2.0d0 )**(gama/(gama-1.0d0))
    enddo
enddo
!
!
do i = 1,points_x
    do j = 1,points_y
       P_total(i,j) = (p_out_trans(i,j)/(gama - 1.0d0)) + Q1(i,j)*(u_out_trans(i,j)**2.0d0+v_out_trans(i,j)**2.0d0)/2.0d0
    enddo
enddo
!
! gravar arquivo para o Tecplot.
!
FileTag = convert(n_save)
!
!
open(UNIT=3,FILE =trim(FileTag)//trim(fname))
    write(3,*)' TITLE="malha do projeto" '
    write(3,*)' VARIABLES = "x", "y" ,"u", "v", "rho", "p", "Mach", "P_estag", "P_Total" '
    write(3,*)' ZONE I=',points_x,' J=',points_y,',DATAPACKING=POINT'
    !
    do j = 1,points_y
       do i = 1,points_x
          WRITE(3,'(7es11.3e2)') meshx(i,j),meshy(i,j),u_out_trans(i,j),v_out_trans(i,j),Q1(i,j),p_out_trans(i,j),M_final(i,j), P_estag(i,j), P_total(i,j)
       enddo
    enddo
    !
    !
close(3)
!
!
!
deallocate(u_out_trans,v_out_trans)
deallocate(p_out_trans,M_final,P_total,P_estag)
!
!
end subroutine output
!
!
!
subroutine output_final
use vars
implicit none
!
!
real(8),dimension(:,:), allocatable         :: u_out_trans, v_out_trans
real(8),dimension(:,:), allocatable         :: p_out_trans, M_final, P_total, P_estag
character nome*10
!
!
!
allocate(u_out_trans(points_x,points_y),v_out_trans(points_x,points_y))
allocate(p_out_trans(points_x,points_y),M_final(points_x,points_y))
allocate(P_total(points_x,points_y),P_estag(points_x,points_y))
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
        p_out_trans(i,j) = (gama -1.0d0)*(Q4(i,j) - Q1(i,j)*(u_out_trans(i,j)**2.0d0+v_out_trans(i,j)**2.0d0)/2.0d0)
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
do i = 1,points_x
    do j = 1,points_y
       P_estag(i,j) = p_out_trans(i,j)*(1.0d0 + ((gama -1.0d0)/2.0d0)*M_final(i,j)**2.0d0 )**(gama/(gama-1.0d0))
    enddo
enddo
!
!
do i = 1,points_x
    do j = 1,points_y
       P_total(i,j) = (p_out_trans(i,j)/(gama - 1.0d0)) + Q1(i,j)*(u_out_trans(i,j)**2.0d0+v_out_trans(i,j)**2.0d0)/2.0d0
    enddo
enddo
!
! gravar arquivo para o Tecplot.
!
open(2,file = 'final.dat')
    !
    write(2,*)' TITLE="malha do projeto" '
    write(2,*)' VARIABLES = "x", "y" ,"u", "v", "rho", "p", "Mach", "P_estag", "P_Total" '
    write(2,*)' ZONE I=',points_x,' J=',points_y,',DATAPACKING=POINT'
    !
    do j = 1,points_y
       do i = 1,points_x
          WRITE(2,'(7es11.3e2)') meshx(i,j),meshy(i,j),u_out_trans(i,j),v_out_trans(i,j),Q1(i,j),p_out_trans(i,j),M_final(i,j), P_estag(i,j), P_total(i,j)
       enddo
    enddo
    !
    !
close(2)
!
!
!
deallocate(u_out_trans,v_out_trans)
deallocate(p_out_trans,M_final,P_total,P_estag)
!
!
end subroutine output_final

