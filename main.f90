PROGRAM verlet_algorithm
  use verlet_module
  IMPLICIT none

  character(len=50) :: file_name
  integer(4) :: n_at, n_traj, ns, i, j, k
  real(8) :: t, dt, kinetic

  character(len=2), dimension(:), allocatable :: atom_list
  real(8), dimension(:,:,:), allocatable :: coord, grad, vel, m
  real(8), dimension(:,:), allocatable :: dummy_cord, mass
  real(8), dimension(:), allocatable :: dummy_mass
  

  ! Initial variables
  file_name="h2o.xyz"
  n_traj=1
  ns=10000
  dt=10.0 !Used to be 10 !a.u
  t=0

  ! Read the input file with intial coordinates and masses 
  call read_input(file_name, atom_list, dummy_cord, dummy_mass, n_at)

  ! Create ensemble arrays
  call create_ensemble_arrays(n_at, n_traj, dummy_mass, dummy_cord, coord, grad, vel, mass)

  ! Convert units
  mass = mass * 1822.888486 !Mass to a.u (electron rest mass)
  coord = coord / 0.52917721092 !Angstrom to Bohr

  ! Obtain the first gradient values
  do i=1,n_traj
     call first_gradient(atom_list, n_at, i, coord(:,:,i), grad(:,:,i))
  end do


  
  !---------------------------------------------------------------------------
  !                               RUN THE DYNAMICS                           !
  !---------------------------------------------------------------------------
  do i=1,ns
     do j=1,n_traj
        
        ! Propagate trajectory j one step
        call MD_step(atom_list, n_at, dt, mass, j, coord(:,:,j), vel(:,:,j), grad(:,:,j))

        ! Write the coordinates after the step was perfomed to a file in order to create the movie
        call trajectory_movie(atom_list, n_at, j, coord(:,:,j), vel(:,:,j))
        
     end do
     
     t = t + dt
     
     kinetic = 0
     do k=1,n_at
        kinetic = kinetic + 0.5*sum(vel(k,:,1)**2)/dummy_mass(k)
     end do
     write(6,*) t,kinetic
      
  end do
  

  
  STOP
END PROGRAM verlet_algorithm
  
