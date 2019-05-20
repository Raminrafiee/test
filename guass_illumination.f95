program guass_ill

  ! use guass_modul
  implicit none
  integer,parameter :: n=3
  integer :: i,j,k
  integer *8 :: temp
  real *8 :: a21
  real *8, dimension(1) :: t3
  real *8,dimension(n,n) ::a
  real *8,dimension(n,n+1) ::aib
  real *8,dimension(n,1) ::b,x

  ! open (1,file='data,txt')
  a=reshape((/2.0_8,-3.0_8,5.0_8,6.0_8,7.0_8,-1.0_8,2.0_8,4.0_8,0.0_8/),(/n,n/))
  b=reshape((/-1.0_8,2.0_8,0.0_8/),(/n,1/))
  print *, "The  matrix A"
  call print_matrics(a)
  print *, "The matrix B"
  call print_matrics(b)


  call augmented_matrix(a,b,n,aib)
  print *, "The augmated matrix"
  call print_matrics(aib)
! -a21/a11 the factor that makes all parameter to zero
do i=1,n-1
  a21=-aib(i+1,1)/aib(1,1)
  do j=1,n+1
    ! print *,j,a21
    aib(i+1,j)=(a21)*a(1,j)+aib(i+1,j)
  end do
end do

print *, "The New matrix"
call print_matrics(aib)
! -a21/a11 the factor that makes all parameter to zero
do i=n,1,-1
t3=matmul(aib(i,1:n),x)
x(i,1)=aib(i,n+1)-t3(1)/aib(i,i)

end do

print *,"Solution"
call print_matrics(x)

contains


        subroutine augmented_matrix (a,b,n,aib)
          implicit none
          integer ,intent(in) :: n
          real *8,dimension(n,n),intent(in) ::a
          real *8,dimension(n),intent(in) ::b
          real *8,dimension(n,n+1),intent(out) ::aib

          aib(1:n,1:n)=a
          aib(1:n,n+1)=b
        end subroutine augmented_matrix

        subroutine print_matrics(a)
          implicit none
          integer :: row,colum
          integer :: i,j
          real(kind=8), dimension(:,:),intent(in) :: a

          row = size(a,dim=1)
          colum = size(a,dim=2)
              5 format(f8.3,t3)
              do i=1,row
                write (*,'("|")',advance='no')
                  do j=1,colum
                      write (*,fmt=5,advance='no') a(i,j)
                  end do
                write (*,'("|")',advance='no')
                write (*,*)
              end do
        end subroutine print_matrics




end program guass_ill
