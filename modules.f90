MODULE verlet_module
  IMPLICIT none

CONTAINS
  SUBROUTINE read_input(file_name, atom_list, coord, mass, n_at)
    !-----------------------------------------------------------
    ! This subroutine reads input files with the following
    ! structure (example for water):
    !
    ! 3
    ! O 0.0000  0.0000  0.0000 16.0000
    ! H 0.0000  0.0000  1.1000 1.000
    ! H 0.6953  0.0000 -0.4345 1.000
    !
    !----------------------------------------------------------
    IMPLICIT NONE
    
    ! Input 
    character(len=50), intent(in) :: file_name ! Input file name
    
    ! Intermediate vars
    integer(4) :: i, j
    
    ! Output
    character(len=2), dimension(:), allocatable, intent(out) :: atom_list
    real(8), dimension(:,:), allocatable, intent(out) :: coord
    real(8), dimension(:), allocatable, intent(out) :: mass
    integer(4), intent(out) :: n_at

    ! Open the input file in unit 20
    open(20,file=file_name)
    read(20,*) n_at

    ! Allocate the arrays
    allocate(atom_list(n_at), coord(n_at,3), mass(n_at))

    ! Read the file
    do i=1,n_at
       read(20,*) atom_list(i), (coord(i,j),j=1,3), mass(i)
    end do

    close(20)
  END SUBROUTINE read_input

  SUBROUTINE create_ensemble_arrays(n_at, n_traj, dummy_mass, dummy_coord, coord, grad, vel, mass)
    !-----------------------------------------------------------
    ! This subroutine transforms data for one trajectory
    ! in arrays that can be used for N trajectories
    !----------------------------------------------------------

    IMPLICIT NONE
    
    ! INPUT
    integer(4), intent(in) :: n_at, n_traj
    real(8), dimension(:,:), intent(in) :: dummy_coord
    real(8), dimension(:), intent(in) :: dummy_mass

    ! INTERMEDIATE
    integer(4) :: i

    ! OUTPUT
    real(8), dimension(:,:,:), allocatable, intent(out) :: coord, grad, vel
    real(8), dimension(:,:), allocatable, intent(out) :: mass


    ! Allocate ensemble arrays
    allocate(coord(n_at,3,n_traj),grad(n_at,3,n_traj),vel(n_at,3,n_traj),mass(n_at,n_at))
    

    ! Create coordinate ensemble array
    do i=1,n_traj
       coord(:,:,i) = dummy_coord(:,:)
    end do

    ! Create N_at x N_at mass array
    do i=1,n_at
       mass(i,:) = dummy_mass(i)
    end do
    

    ! For now, grad and vel ensemble array are set to zero
    grad = 0
    vel = 0
  END SUBROUTINE create_ensemble_arrays


  SUBROUTINE first_gradient(atom_list, n_at, traj, coord, grad)
    !-----------------------------------------------------------
    ! Read the first gradient that will make the nucleii move 
    ! in the field of the electrons. Basically it corresponds
    ! to the dH/dR_ij where i makes reference to the nucleus i
    ! and j to the direction x y or z.
    !----------------------------------------------------------

    IMPLICIT NONE
    
    ! INPUT
    character(len=2), dimension(:), intent(in) :: atom_list
    integer(4), intent(in) :: n_at, traj
    real(8), dimension(:,:), intent(in) :: coord

    ! INTERMEDIATE
    character(len=50) :: run_command
    
    ! OUTPUT
    real(8), dimension(:,:), intent(out) :: grad

    ! Write the geometry of thde trajectory Nj to a file geom_Nj.xyz
    call write_geometry(traj, n_at, atom_list, coord)

    ! Call the system to calculate the gradient
    write(run_command, '(A9,I2.2,A5,I2.2,A5,I2.2)') './run.sh ', traj, " geom", traj,".xyz ", n_at*3
    call system(run_command)
    
    ! Read the gradient
    call read_gradient(n_at, traj, grad)
    
  END SUBROUTINE first_gradient


  SUBROUTINE write_geometry(traj, n_at, atom_list, coord)
    !-----------------------------------------------------------
    ! Writes the current geometry of the trajectory with number
    ! 'traj' to a file named geom'traj'.xyz
    !----------------------------------------------------------

    IMPLICIT NONE
    
    ! INPUT
    character(len=2), dimension(:), intent(in) :: atom_list
    integer(4), intent(in) :: n_at, traj
    real(8), dimension(:,:), intent(in) :: coord
    
    
    ! INTERMEDIATE
    integer(4) :: k, m, io_unit
    character(len=50) :: file_name

    io_unit = traj*2

    write(file_name,'(A4,I2.2,A4)') 'geom', traj, '.xyz'
    open(io_unit,file=file_name)

    write(io_unit, '(I3)') n_at
    write(io_unit,*) 'Angstrom'

    do k = 1, n_at
       write(io_unit,*) atom_list(k), (coord(k,m)*0.52917721092,m=1,3) ! Bohr to Angstrom
    end do

    close(io_unit)
    
        
  END SUBROUTINE write_geometry


  SUBROUTINE read_gradient(n_at, traj, grad)
    !-----------------------------------------------------------
    ! Reads the energy of 'traj' from a file named 'traj'.engrad
    ! Assumes that gradiends are in atomic units (Ha/Bohr)
    !-----------------------------------------------------------

    IMPLICIT NONE
    
    ! INPUT
    integer(4), intent(in) :: n_at, traj


    ! INTERMEDIATE
    integer(4) :: k, m, io_unit
    character(len=50) :: file_name

    ! OUTPUT
    real(8), dimension(:,:), intent(out) :: grad

    io_unit = traj*2

    write(file_name,'(A4,I2.2,A7)') 'traj', traj, '.engrad'
    open(io_unit, file=file_name)

    ! Read the gradient form traj01.engrad file (for example)
    do k=1,n_at
       do m=1,3
          read(io_unit,*) grad(k,m) 
       end do   
    end do

       
  END SUBROUTINE read_gradient
  
  SUBROUTINE MD_step(atom_list, n_at, dt, mass, traj, coord, vel, grad)
    !-----------------------------------------------------------
    ! Propagates the nucleii of a certrain trajectory
    ! using the gradients obtained from
    ! electronic calculations. Therefore, the nucleii will move
    ! on the field of the electrons.
    ! The velocity-verlet algorithm is used to do the numerical
    ! integration.
    !-----------------------------------------------------------

    IMPLICIT NONE
    
    ! INPUT
    integer(4), intent(in) :: n_at, traj
    real(8), intent(in) :: dt
    real(8), dimension(:,:), intent(in) :: mass
    character(len=2), dimension(:), intent(in) :: atom_list

    ! INTERMEDIATE
    character(len=50) :: run_command

    ! OUTPUT
    real(8), dimension(:,:), intent(out) :: coord, vel, grad
    
    ! Calculate Vel(n+1/2) and Coord(n+1)
    vel = vel + 0.5 * dt * (-grad/mass) 
    coord = coord + vel * dt

    ! Write the geometry of the current trajectory Nj to a file geom_Nj.xyz
    call write_geometry(traj, n_at, atom_list, coord)
    
    ! Call the system to calculate the new gradient Grad(n+1)
    write(run_command, '(A9,I2.2,A5,I2.2,A5,I2.2)') './run.sh ', traj, " geom", traj,".xyz ", n_at*3
    call system(run_command)
    
    ! Read the new gradients Grad(n+1)
    call read_gradient(n_at, traj, grad)

    ! Calculate Vel(n+1)
    vel = vel + 0.5 * dt * (-grad/mass)
    
  END SUBROUTINE MD_step
  

  
  SUBROUTINE trajectory_movie(atom_list, n_at, traj, coord, vel)
    !-----------------------------------------------------------
    ! Write the trajectory coordinates to a file in order to
    ! see the trajectory in molde (for example).
    !-----------------------------------------------------------

    IMPLICIT NONE
    
    ! INPUT
    integer(4), intent(in) :: n_at, traj
    character(len=2), dimension(:), intent(in) :: atom_list
    real(8), dimension(:,:), intent(out) :: coord, vel

    ! INTERMEDIATE
    character(len=50) :: file_name
    integer(4) :: k, m, io_unit
   

    io_unit = traj*2
    
    write(file_name,'(A4,I2.2,A4)') 'traj', traj, '.xyz'
    
    open(io_unit,position="append",file=file_name)
    write(io_unit, '(I3)') n_at
    write(io_unit,*) 'Angstrom'

    do k = 1, n_at
       write(io_unit,*) atom_list(k), (coord(k,m)*0.529177249,m=1,3) ! Bohr to Angstrom
    end do

    close(io_unit)

    
  END SUBROUTINE trajectory_movie
  

  
END MODULE verlet_module
