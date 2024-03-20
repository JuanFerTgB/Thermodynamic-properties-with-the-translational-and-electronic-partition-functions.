program  Gibbs_free_energy_entropy_internal_energy   

    implicit none
    double precision, parameter :: e = 2.7182818284590452353602874719777d0                                          ! Euler's number
    double precision, parameter :: pi = 3.14159265358979323846d0                                                    ! Pi
    double precision :: kb , h ,P ,V ,n, R,  T ,  a_v , U , R2 ,geg  ,qtus,entropi_eus,entropius,LG_tus,LG_eus,entropi_teus,GF_teus    ! Variables of the general problem.
    real :: m_us
    real(8) :: m(9),  qt(9), entropi(9), entropi_e(9), ge(9), LG_t(9) , LG_e(9) , entropi_te(9) , GF_te(9)          ! Tools for the do loop
    integer :: i                                                                                                    ! Variable i
    open(100, file='values.text' , status='replace'  )
    ! Constants
     kb  = 1.38064852e-23                   ! Boltzmann constant (J/K)
     a_v = 6.022e23                         ! Avogadros Number
     h   = 6.62607004e-34                   ! Planck constant (J s)
     R2 = 8.3145D0                          ! Ideal gas constant
    ! Masses in kg of the atoms
     m(1) = 4.003 *(1.0/1000.0) *(1/a_v)    ! Helium mass (kg)
     m(2) = 20.18 *(1.0/1000.0) *(1/a_v)    ! Neon mass (kg)
     m(3) = 26.98*(1.0/1000.0) *(1/a_v)     ! Aluminum mas (kg)
     m(4) = 137.34 *(1.0/1000.0) *(1/a_v)   ! Barium mass (kg)
     m(5) = 39.94*(1.0/1000.0) *(1/a_v)     ! Argon mass (kg)
     m(6) = 9.01*(1.0/1000.0) *(1/a_v)      ! Beryllium mass (kg)
     m(7) = 159.80*(1.0/1000.0) *(1/a_v)    ! Bromine mass (kg)
     m(8) = 12.01 *(1.0/1000.0) *(1/a_v)    ! Carbon mass (kg)
     m(9) = 1.008  *(1.0/1000.0) *(1/a_v)   ! Hydrogen mass (kg)
     T = 298.15d0          ! Temperature (K)
     P= 1.0                ! Bar pressure 
     n= a_v                 ! N_A
     R= 0.08205746d0        ! Ideal gas constant
     
    ! To solve for the molar volume
    V = ((R*T)/n*P)/1000.0
    
    ! Calculation of the translational partition function
    do i = 1, 9
    qt(i) = ((2.0 * pi * m(i) * kb * T) / h**2)**(3.0/2.0) * V
    end do
    

    ! The electronic partition function for atoms is almost equal to the degeneracies
    ! Degeneracies for the calculation of entropy 
    ge(1) = 1 ! Helium 
    ge(2) = 1 ! Neon 
    ge(3) = 2 ! Aluminum 
    ge(4) = 1 ! Barium 
    ge(5) = 1 ! Argon 
    ge(6) = 1 ! Beryllium 
    ge(7) = 4 ! Bromine 
    ge(8) = 1 ! Carbon 
    ge(9) = 1 ! Hydrogen 

    call entropy_translational                                  ! Calculation of the entropy using the translational partition function
    call entropy_electronic                                     ! Calculation of the entropy using the electronic partition function
    call Gibbs_free_energy_translational                        ! Calculation of the Gibbs free energy using the translational partition function
    call Gibbs_free_energy_electronic                           ! Calculation of the Gibbs free energy using the electronic partition function
    call entropy_translational_and_electronic                   ! Calculation of entropy using both partition function
    call Gibbs_free_energy_translational_and_electronic         ! Calculation of Gibbs free energy using both partition function
    call Internal_energy                                        ! Calculation of internal energy 
    call Electronic_configuration                               ! Selection of the orbital and degeneracy.
    call User_input_calculation                                 ! User input
    
contains

    subroutine entropy_translational
       
        ! Calculation of entropy using the translational partition function
        do i = 1, 9
        entropi(i) = kb * log(qt(i) * exp(3.0/2.0)) * a_v
        end do
        
        ! Print the results of the entropy calculated using the translational partition function
        do i = 1, 9
            select case (i)
                case (1)
                    write(100, '(A,F6.2,A)') "Translational entropy of Helium:", entropi(i), "J/mol*K"
                case (2)
                    write(100, '(A,F6.2,A)') "Translational entropy of Neon:", entropi(i), "J/mol*K"
                case (3)
                    write(100, '(A,F6.2,A)') "Translational entropy of Aluminum:", entropi(i), "J/mol*K"
                case (4)
                    write(100, '(A,F6.2,A)') "Translational entropy of Barium:", entropi(i), "J/mol*K"
                case (5)
                    write(100, '(A,F6.2,A)') "Translational entropy of Argon:", entropi(i), "J/mol*K"
                case (6)
                    write(100, '(A,F6.2,A)') "Translational entropy of Beryllium:", entropi(i), "J/mol*K"
                case (7)
                    write(100, '(A,F6.2,A)') "Translational entropy of Bromine:", entropi(i), "J/mol*K"
                case (8)
                    write(100, '(A,F6.2,A)') "Translational entropy of Carbon:", entropi(i), "J/mol*K"
                case (9)
                    write(100, '(A,F6.2,A)') "Translational entropy of Hydrogen:", entropi(i), "J/mol*K"
            end select
        end do
           
    end subroutine entropy_translational

    subroutine  entropy_electronic

        ! Calculation of entropy using the electronic partition function
        do i = 1, 9
        entropi_e(i) = kb * log(ge(i) ) * a_v
        end do
        
        ! Print the results of the entropy calculated using the electronic partition function
        do i = 1, 9
            select case (i)
                case (1)
                    write(100, '(A,F6.2,A)') "Electronic entropy of Helium:", entropi_e(i), "J/mol*K"
                case (2)
                    write(100, '(A,F6.2,A)') "Electronic entropy of Neon:", entropi_e(i), "J/mol*K"
                case (3)
                    write(100, '(A,F6.2,A)') "Electronic entropy of Aluminum:", entropi_e(i), "J/mol*K"
                case (4)
                    write(100, '(A,F6.2,A)') "Electronic entropy of Barium:", entropi_e(i), "J/mol*K"
                case (5)
                    write(100, '(A,F6.2,A)') "Electronic entropy of Argon:", entropi_e(i), "J/mol*K"
                case (6)
                    write(100, '(A,F6.2,A)') "Electronic entropy of Beryllium:", entropi_e(i), "J/mol*K"
                case (7)
                    write(100, '(A,F6.2,A)') "Electronic entropy of Bromine:", entropi_e(i), "J/mol*K"
                case (8)
                    write(100, '(A,F6.2,A)') "Electronic entropy of Carbon:", entropi_e(i), "J/mol*K"
                case (9)
                    write(100, '(A,F6.2,A)') "Electronic entropy of Hydrogen:", entropi_e(i), "J/mol*K"
            end select
        end do

    end subroutine entropy_electronic

    subroutine Gibbs_free_energy_translational

        ! Calculation of Gibbs free energy using the translational partition function
        do i = 1,9
            LG_t(i) = -kb*T*log(qt(i))*a_v + kb*T*a_v
        end do
        
        ! Print the results of the Gibbs free energy calculated using the translational partition function
        do i = 1, 9
            select case (i)
                case (1)
                    write(100, '(A,F10.2,A)') "Translational Gibbs free energy of Helium:", LG_t(i), "J"
                case (2)
                    write(100, '(A,F10.2,A)') "Translational Gibbs free energy of Neon:", LG_t(i), "J"
                case (3)
                    write(100, '(A,F10.2,A)') "Translational Gibbs free energy of Aluminum:", LG_t(i), "J"
                case (4)
                    write(100, '(A,F10.2,A)') "Translational Gibbs free energy of Barium:", LG_t(i), "J"
                case (5)
                    write(100, '(A,F10.2,A)') "Translational Gibbs free energy of Argon:", LG_t(i), "J"
                case (6)
                    write(100, '(A,F10.2,A)') "Translational Gibbs free energy of Beryllium:", LG_t(i), "J"
                case (7)
                    write(100, '(A,F10.2,A)') "Translational Gibbs free energy of Bromine:", LG_t(i), "J"
                case (8)
                    write(100, '(A,F10.2,A)') "Translational Gibbs free energy of Carbon:", LG_t(i), "J"
                case (9)
                    write(100, '(A,F10.2,A)') "Translational Gibbs free energy of Hydrogen:", LG_t(i), "J"
            end select
        end do

    end subroutine Gibbs_free_energy_translational

    subroutine Gibbs_free_energy_electronic
        
        ! Calculation of Gibbs free energy using the electronic partition function
        do i = 1,9
            LG_e(i) = -kb*T*log(ge(i))*a_v + (kb*T*a_v)
        end do
        
        ! Print the results of the Gibbs free energy calculated using the electronic partition function
        do i = 1, 9
            select case (i)
                case (1)
                    write(100,'(A,F10.2,A)') "Electronic Gibbs free energy of Helium:", LG_e(i), "J"
                case (2)
                    write(100, '(A,F10.2,A)') "Electronic Gibbs free energy of Neon:", LG_e(i), "J"
                case (3)
                    write(100, '(A,F10.2,A)') "Electronic Gibbs free energy of Aluminum:", LG_e(i), "J"
                case (4)
                    write(100, '(A,F10.2,A)') "Electronic Gibbs free energy of Barium:", LG_e(i), "J"
                case (5)
                    write(100, '(A,F10.2,A)') "Electronic Gibbs free energy of Argon:", LG_e(i), "J"
                case (6)
                    write(100, '(A,F10.2,A)') "Electronic Gibbs free energy of Beryllium:", LG_e(i), "J"
                case (7)
                    write(100, '(A,F10.2,A)') "Electronic Gibbs free energy of Bromine:", LG_e(i), "J"
                case (8)
                    write(100, '(A,F10.2,A)') "Electronic Gibbs free energy of Carbon:", LG_e(i), "J"
                case (9)
                    write(100, '(A,F10.2,A)') "Electronic Gibbs free energy of Hydrogen:", LG_e(i), "J"
            end select
        end do

    end subroutine Gibbs_free_energy_electronic

    subroutine entropy_translational_and_electronic
    
        ! Calculation of entropy using both the electronic and translational partition functions

            do i =1,9
                !entropi_te(i) = kb * log(qt(i) * exp(3.0/2.0)*ge(i)) * a_v
                entropi_te(i)= (5.0/2.0)*R2 + R2*log(   qt(i) * ge(i)   )
            end do

        ! Print the results of entropy using both the electronic and translational partition functions
            do i = 1, 9
                select case (i)
                    case (1)
                        write(100, '(A,F6.2,A)') "Entropy using both partition functions of Helium:", entropi_te(i), "J/mol*K"
                    case (2)
                        write(100, '(A,F6.2,A)') "Entropy using both partition functions of Neon:", entropi_te(i), "J/mol*K"
                    case (3)
                        write(100, '(A,F6.2,A)') "Entropy using both partition functions of Aluminum:", entropi_te(i), "J/mol*K"
                    case (4)
                        write(100, '(A,F6.2,A)') "Entropy using both partition functions of Barium:", entropi_te(i), "J/mol*K"
                    case (5)
                        write(100, '(A,F6.2,A)') "Entropy using both partition functions of Argon:", entropi_te(i), "J/mol*K"
                    case (6)
                        write(100, '(A,F6.2,A)') "Entropy using both partition functions of Beryllium:", entropi_te(i), "J/mol*K"
                    case (7)
                        write(100, '(A,F6.2,A)') "Entropy using both partition functions of Bromine:", entropi_te(i), "J/mol*K"
                    case (8)
                        write(100, '(A,F6.2,A)') "Entropy using both partition functions of Carbon:", entropi_te(i), "J/mol*K"
                    case (9)
                        write(100, '(A,F6.2,A)') "Entropy using both partition functions of Hydrogen:", entropi_te(i), "J/mol*K"
                end select
            end do

    end subroutine entropy_translational_and_electronic

    subroutine Gibbs_free_energy_translational_and_electronic     
    
        ! Calculation of Gibbs free energy using using both the electronic and translational partition functions
            do i = 1,9
                GF_te(i) = -kb*a_v*T*log(ge(i)*qt(i)/a_v)+kb*T*a_v
            end do

        ! Print the results of the Gibbs free energy calculated using both the electronic and translational partition functions
            do i = 1, 9
                select case (i)
                    case (1)
                        write(100, '(A,F10.2,A)') "Gibbs free energy  using both partition functions of Helium:",GF_te(i), "J"
                    case (2)
                        write(100, '(A,F10.2,A)') "Gibbs free energy  using both partition functions of Neon:", GF_te(i), "J"
                    case (3)
                        write(100, '(A,F10.2,A)') "Gibbs free energy  using both partition functions of Aluminum:", GF_te(i), "J"
                    case (4)
                        write(100, '(A,F10.2,A)') "Gibbs free energy  using both partition functions of Barium:", GF_te(i), "J"
                    case (5)
                        write(100, '(A,F10.2,A)') "Gibbs free energy  using both partition functions of Argon:", GF_te(i), "J"
                    case (6)
                        write(100, '(A,F10.2,A)') "Gibbs free energy  using both partition functions of Beryllium:", GF_te(i), "J"
                    case (7)
                        write(100, '(A,F10.2,A)') "Gibbs free energy  using both partition functions of Bromine:", GF_te(i), "J"
                    case (8)
                        write(100, '(A,F10.2,A)') "Gibbs free energy  using both partition functions of Carbon:", GF_te(i), "J"
                    case (9)
                        write(100, '(A,F10.2,A)') "Gibbs free energy  using both partition functions of Hydrogen:", GF_te(i), "J"
                end select
            end do

        
    end subroutine Gibbs_free_energy_translational_and_electronic     

    subroutine Electronic_configuration

        character(len=100) :: electron_config
        character(len=10) :: last_term
        integer ::  len_config, orbital_start
    
        ! Read the electronic configuration.
    
        print*, "Enter the electronic configuration example 1s^2 2p^2:"
        read(*, '(A)') electron_config
    
        ! Find the length of the string.
    
        len_config = len_trim(electron_config)
    
        ! Find the term referring to the orbital.
    
        orbital_start = len_config
        do i = len_config, 1, -1
            if (electron_config(i:i) == ' ') exit
            if (electron_config(i:i) == '^') orbital_start = i - 1
        end do
    
        ! Extract the last term (orbital).
    
        last_term = electron_config(orbital_start:len_config)
    
        !Print the last term (orbital).
    
        print*, "The orbital is:", last_term

        
        ! The degeneracy in this configuration is 1
        if (last_term == 's^2') then
            geg=1
        elseif (last_term == 'p^6') then
            geg= 1
        elseif (last_term == 'd^10') then
            geg=1
        elseif (last_term == 'p^2') then
            geg=1
        elseif (last_term == 'd^4') then
            geg=1

        ! The degeneracy in this configuration is 2
        elseif (last_term == 'p^1') then
            geg= 2
        
        ! The degeneracy in this configuration is 4
        elseif (last_term == 'p^5') then
            geg=4
        elseif (last_term == 'd^1') then
            geg=4
        elseif (last_term == 'd^3') then
            geg=4

        ! The degeneracy in this configuration is 5
        elseif (last_term == 'p^4') then
            geg= 5
        elseif (last_term == 'd^2') then
            geg=5

        ! The degeneracy in this configuration is 6
        elseif (last_term == 'p^3') then
            geg=6
        elseif (last_term == 'd^9') then
            geg=6

        ! The degeneracy in this configuration is 9
        elseif (last_term == 'd^8') then
            geg=9
        elseif (last_term == 'd^6') then
            geg=9

        ! The degeneracy in this configuration is 10
        elseif (last_term == 'd^7') then
            geg= 10
        elseif (last_term == 'd^5') then
            geg=10
        endif
    
    end subroutine Electronic_configuration

    subroutine Internal_energy 
        ! Calculation of the internal energy, which is constant for all atoms
        U = 3.0/2.0 * a_v *kb *T
        ! Print the results of theinternal energy
            write(100, '(A,F10.2,A)') "Internal energy using both partition functions", U, "J"       
    end subroutine Internal_energy

    subroutine User_input_calculation
        ! Request for the atom's mass from the user
        write(*,*) 'Enter the value of the atoms mass'
        read(*,*) m_us

        ! Translational function of the user atom
        qtus = ((2.0 * pi *m_us * kb * T) / h**2)**(3.0/2.0) * V  
        write(100,'(A,F10.2,A)') 'Translational function of the user atom',qtus, ''

        ! Translational entropy of the users atom
        entropius = kb * log(qtus * exp(3.0/2.0)) * a_v
        write(100,'(A,F10.2,A)') 'Translational entropy of the users atom',entropius, 'J/mol*K'

        ! Electronic entropy of the users atom
        entropi_eus = kb * log(geg ) * a_v
        write(100,'(A,F10.2,A)') 'Electronic entropy of the users atom',entropi_eus, 'J/mol*K'

        ! Gibbs free energy of translation of the users atom
        LG_tus = -kb*T*log(qtus)*a_v + kb*T*a_v
        write(100,'(A,F10.2,A)') 'Gibbs free energy of translation of the users atom',LG_tus, 'J'

        ! Gibbs free energy of electronic of the users atom
        LG_eus = -kb*T*log(geg)*a_v + (kb*T*a_v)
        write(100,'(A,F10.2,A)') 'Gibbs free energy of electronic of the users atom',LG_eus, 'J'
        
        ! Entropy using both Partition functions of the users atom
        entropi_teus= (5.0/2.0)*R2 + R2*log(   qtus * geg   )
        write(100,'(A,F10.2,A)') 'Entropy using both Partition functions of the users atom',entropi_teus, 'J/mol*K'

        ! Gibbs free energy using both Partition functions of the users atom
        GF_teus = -kb*a_v*T*log(geg*qtus/a_v)+kb*T*a_v
        write(100,'(A,F10.2,A)') 'Gibbs free energy using both Partition functions of the users atom ',GF_teus, 'J'
    end subroutine User_input_calculation

end program Gibbs_free_energy_entropy_internal_energy   
