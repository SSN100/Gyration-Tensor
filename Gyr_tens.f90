Program Gyr_tensor
    implicit none
    !Declaration
    integer, parameter :: N=16720
    integer :: i, n_frames=1250, j
    real :: coord(3, N), coord_M(3, N), mass(N), CM(3), gyr_ten(3, 3)
    character :: atom_name(N)
    
    open(10, file="npt_prod_skip_4_not_sc_only_coord.dat")
    open(11, file="atom_name_new.dat")
    open(12, file="Gyration_tensor_2.dat")
    open(13, file="Center_of_Mass.dat")
    
    !Read a file containing atom_name which corresponding coordinate can be read later.
    read(11, *) atom_name
    
    
    do j=1, n_frames
    
        coord=0.0
        do i=1, N
            read(10, *) coord(:, i) !Read a file which contains only coordinate.
            select case (atom_name(i))
		!Assign mass to every atom
                case ('S')
                    mass(i)=32.065
                case ('O')
                    mass(i)=15.999
                case ('N')
                    mass(i)=14.007
                case ('C')
                    mass(i)=12.01
                case ('H')
                    mass(i)=1.008
            end select
            coord_M(:, i)=coord(:, i)*mass(i)
        end do
	
	!The coordinate file has 3 unnecessary lines between every frames and hence skip them
        do i=1, 3
            read(10, *)
        end do
    
    
        CM=sum(coord_M, dim=2)/sum(mass) ! Getting Center of Mass
!         write(13, "(i0)") j
        write(13, "(f15.5, 3f15.5)") real(j)*4.0*10.0, CM
        

	!Getting gyration tensor.
        gyr_ten=0.0
        do i=1, N
            gyr_ten(1, 1)=gyr_ten(1, 1)+(coord(1, i)-CM(1))**2
            gyr_ten(2, 1)=gyr_ten(2, 1)+(coord(1, i)-CM(1))*(coord(2, i)-CM(2))
            gyr_ten(3, 1)=gyr_ten(3, 1)+(coord(1, i)-CM(1))*(coord(3, i)-CM(3))
            gyr_ten(1, 2)=gyr_ten(1, 2)+(coord(2, i)-CM(2))*(coord(1, i)-CM(1))
            gyr_ten(2, 2)=gyr_ten(2, 2)+(coord(2, i)-CM(2))**2
            gyr_ten(3, 2)=gyr_ten(3, 2)+(coord(2, i)-CM(2))*(coord(3, i)-CM(3))
            gyr_ten(1, 3)=gyr_ten(1, 3)+(coord(3, i)-CM(3))*(coord(1, i)-CM(1))
            gyr_ten(2, 3)=gyr_ten(2, 3)+(coord(3, i)-CM(3))*(coord(2, i)-CM(2))
            gyr_ten(3, 3)=gyr_ten(3, 3)+(coord(3, i)-CM(3))**2
        end do
        gyr_ten=gyr_ten/N
!         write(12, "(i0)") j
        !write(12, "(f15.5)") real(j)*4.0*10.0
        write(12, "(3f15.5)") gyr_ten
        
    end do
    
End Program
