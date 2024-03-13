program planets
    use processes1
    implicit none
    real*8, parameter :: X_sample = 4320d0, Y_sample = 4320d0
    real*8, parameter :: d = 0.27d0, ErgToEv = 6.241509651616162d+11, F0 = 2.07d-7 
    real*8, parameter :: lambda_sgs = 1.2d-5, lamb_sgs_inv = 8.333333d+4 !эрстеды
    real*8:: A, B, energy, K0_const, Def_const, H0, B_const, B_sum
    real :: start_time, finish_time
    integer:: n, cnt
    K0_const = ((F0**2)*(lamb_sgs_inv**3))/(8d0*((3.14d0)**2)) 
    Def_const = -0.8d0*K0_const
    B_const = (F0*(lamb_sgs_inv**2))/(2*3.14d0)
    cnt = 0
    A=0.0
    B = 100000.0
    n=10000
    print*, K0_const
    call cpu_time(start_time)
    call init_config(energy)
    H0 = 1000d0
    call DE10(A, B, n)
    call cpu_time(finish_time)
    write(*,'("Run time = ", F10.5)') finish_time-start_time
    contains

    subroutine DE10(A1, B1, N1)
      real*8, intent(in):: A1, B1
      integer, intent(in):: N1
      !real*8, intent(inout):: Y1(:)
      real*8, parameter:: C0(3)=(/0.5, 0.5, 1.0/),   &  
        P(3)=(/0.16666667, 0.33333333, 0.33333333/)
      real*8:: CH(3)
      real*8, allocatable:: Y0(:), YH(:), K(:), Y0_old(:)
      real*8:: X, XH, H
      integer:: i, j, l, M1, i3, z, M_old
        M1=size(Y)
        do i3 = 1, M1-1, 2
          open(11, file='p0.txt')
          write(11,'(F36.16, ",", F36.16)') Y(i3), Y(i3+1) !Y1(i3), Y1(i3+1)
        end do
        if (A1==B1) return
        H=(B1-A1)/float(N1) !Задаем шаг сетки по времени
        do j=1, 3
          CH(j)=H*C0(j)
        end do
        allocate (Y0(M1), YH(M1), K(M1), Y0_old(M1))
        Y0=Y
        Y0_old = Y0
        M_old = M1
        do l=0, N1-1
          X=A1+H*float(l) 
          call DS(H0, Y, K) !K здесь перезаписывается по новому
          do j=1, 3
            XH=X+CH(j)  
            do i=1, M1
              K(i)=H*K(i) 
              YH(i)=Y0(i)+C0(j)*K(i)  !YH тут тоже по новому
              Y(i)=Y(i)+P(j)*K(i)
            end do
            call DS(H0, YH, K)
          end do
          do i=1, M1
            Y(i)=Y(i)+P(1)*H*K(i) 
          end do
          !Полностью периодическая граница
          do i = 1, M1-1, 2
            if (Y(i+1) >= Y_sample) then
              Y(i+1) = Y(i+1) - Y_sample
            endif
            if (Y(i+1) <= 0d0) then   
              Y(i+1) = Y(i+1) + Y_sample 
            endif
            if (Y(i) >= X_sample) then
              Y(i) = Y(i) - X_sample 
            endif
            if (Y(i) <= 0d0) then
              Y(i) = Y(i) + X_sample
            end if
          end do
          Y0=Y
          if (mod(l,10) == 0) then
              do i3 = 1, M1-1, 2
                open(11, file='p'//trim(str(l+10))//'.txt')
                write(11,'(F36.16, ",", F36.16)') Y(i3), Y(i3+1)
              end do
          endif
          if ((l <= N1/2 .AND. l >= 450) .AND. mod(l,450) == 0) then 
            call add_vortices(energy)
            M1 = size(Y)
            H0 = H0 + 100d0
            deallocate (Y0, YH, K)
            allocate (Y0(M1), YH(M1), K(M1))
            do z = 1, M_old, 1
              Y0(z) = Y0_old(z)
            end do
            do z = M_old+1, M1, 1
              Y0(z) = Y(z)
            end do
            deallocate(Y0_old)
            allocate(Y0_old(M1))
            Y0_old = Y0
            M_old = M1
          end if
          if (l > N1/2 .AND. mod(l,450) == 0) then
            call delete_vortices(energy)
            M1 = size(Y)
            H0 = H0 - 100d0
            deallocate (Y0, YH, K)
            allocate (Y0(M1), YH(M1), K(M1))
            do z = 1, M1, 1
              Y0(z) = Y0_old(z)
            end do
            deallocate(Y0_old)
            allocate(Y0_old(M1))
            Y0_old = Y0
            M_old = M1
          end if
          if (l >= 400 .AND. mod(l,400)==0) then
            call B_calc(Y, B_sum)
            B_sum = B_const*B_sum
            open(32, file='B_average.txt')
            write(32,'(F36.16, ",", F36.16)') H0, B_sum
            !print*, B_sum
          endif
        end do
        deallocate (Y0, YH, K)
        return
      end subroutine DE10
      

      subroutine  DS(H, Y2, DY)
        implicit none
        real*8 :: H, Y2(:)
        real*8:: DY(:)
        real*8 :: DY_x, DY_y, X_defs, Y_defs, Meissner_const
        integer :: j2, M2 
        Meissner_const = (F0*H*lamb_sgs_inv)/(4d0*(3.14d0)*dcosh(0.5d0*X_sample*lambda_inv))
        M2 = size(Y2)
        DY_x = 0d0
        DY_y = 0d0
        do j2 = 1, M2-1, 2
          call interaction_vortices(j2, Y2, DY_x, DY_y) !взаимод одного вихря со всеми
          call int_wth_defects(j2, Y2, X_defs, Y_defs) !взаимод одного вихря со всеми дефектами
          DY(j2) = K0_const*DY_x - Meissner_const * dsinh((Y2(j2)-(X_sample*0.5d0))*lambda_inv) + Def_const*X_defs
          DY(j2+1) = K0_const*DY_y + Def_const*Y_defs
        end do
        return
    end subroutine DS


        character(len=20) function str(k)
        !   "Convert an integer to string."
          integer, intent(in) :: k
          write (str, *) k
          str = adjustl(str)
        end function str
end program planets