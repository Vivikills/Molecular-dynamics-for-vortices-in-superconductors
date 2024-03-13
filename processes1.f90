module processes1
    use rand
    implicit none
    real*8, allocatable :: Y(:), YD(:)
    integer :: Pcle_Num, num_rand1, num_rand2
    real*8, parameter :: T = 4.2d0
    real*8, parameter :: x_grid = 2.0d0, lambda_inv = 0.008333333d0
    real*8, parameter :: X_size = 4320d0, Y_size = 4320d0, sigma = 0.6d0, esh = 10d0 !здесь сигма уже в ангстремах
    contains

    subroutine init_config(system_energy)
        implicit none
        integer :: i, j, count, Numb_def, m, j1, i1, t1, count1
        real*8 :: dis, x_part, y_part, yy
        real*8, intent(inout) :: system_energy
        num_rand1 = 4839
        num_rand2 = 29472
        call SRANMAR(num_rand1,  num_rand2)
        Pcle_Num = 900
        Numb_def = 990
        write(*,'("Initial number of particles: ", I3)') Pcle_Num
        allocate(Y(2*Pcle_Num))
        allocate(YD(2*Numb_def))
        count = 0
        system_energy = 0d0
        ! Square distribution vortices
        ! do i = 1, 10, 1
        !     y_part = 300d0*(i-1) + 150d0
        !     count = 0
        !     do j = 1 + (i-1)*6, i*6, 1
        !         count = count + 1
        !         Y(j*2-1) = 300d0*(count-1) + 750d0     
        !         Y(j*2) = y_part
        !         write(*, '("X", I3, " = ", F7.1, ", Y", I3, " = ", F7.1)') j, Y(j*2-1), j, Y(j*2)
        !     enddo
        ! enddo
        ! Random distribution on 120 vortices
        do i = 1, 2*Pcle_Num - 1, 2  
            count = count + 1      
            if (i <= Pcle_Num) x_part = 120d0* RNDM()       
            if (i > Pcle_Num) x_part = 120d0* (RNDM()-1) + X_size
            y_part = Y_size * RNDM()
            Y(i) = x_part         
            Y(i+1) = y_part  
            do j = i-2, 1, -2
                call distance(x_part, y_part, Y(j), Y(j+1), dis)
                if ((x_part == Y(j)  .AND.  y_part == Y(j+1)) .OR.&
                    (dis == 2)) then       
                    print*, "Rechoosing coordinates"              
                    x_part = 120d0 * RNDM()
                    Y(i) = x_part                 
                endif
            end do
            write(*, '("X", I3, " = ", F7.1, ", Y", I3, " = ", F7.1)') i, x_part, i, y_part
        end do
        !Random distribution vortices
        ! do i = 1, 2*Pcle_Num-1, 2  
        !   count = count + 1      
        !   x_part = RNDM()*X_size
        !   y_part = Y_size * RNDM()
        !   Y(i) = x_part         
        !   Y(i+1) = y_part  
        !   do j = i-2, 1, -2
        !       call distance(x_part, y_part, Y(j), Y(j+1), dis)
        !       if ((x_part == Y(j)  .AND.  y_part == Y(j+1)) .OR.&
        !           (dis <= 2)) then       
        !           print*, "Rechoosing coordinates"              
        !           x_part = X_size*RNDM()
        !           Y(i) = x_part                 
        !       endif
        !   end do
        !   write(*, '("X", I3, " = ", F7.1, ", Y", I3, " = ", F7.1)') i, x_part, i, y_part
        ! end do
        !open(14, file='def_coords2.txt')
        !Random distribution defects
        ! do i = 1, 2*Numb_def-1, 2
        !     YD(i) = RNDM()*(X_size-240d0) + 120d0
        !     YD(i+1) = Y_size*RNDM() 
        !     do j1 = i-2, 1, -2
        !       call distance(YD(i), YD(i+1), YD(j1), YD(j1+1), dis1)
        !       !print*, dis1
        !       if (dis1 < 64d0) then 
        !         print*, "Rechoosing coordinates"              
        !         YD(i+1) = Y_size*RNDM()               !неправильный цикл
        !       endif
        !     end do
        !     write(12,'(F36.16, ",", F36.16)') YD(i), YD(i+1)
        !     write(*, '("X", I3, " = ", F7.1, ", Y", I3, " = ", F7.1)') i, YD(i), i, YD(i+1)
        ! enddo
        !square distribution defects
        ! do i = 1, 10, 1
        !     YD(i+1) = 120d0*(i-1) + 60d0 ! снизу вверх заполнение 
        !     count = 0
        !     do j = 1 + (i-1)*10, i*10, 1
        !         count = count + 1
        !         YD(j*2-1) = 120d0*(count-1) + 180d0     
        !         YD(2*j) = YD(i+1)
        !         write(*, '("X", I3, " = ", F7.1, ", Y", I3, " = ", F7.1)') j, YD(j*2-1), j, YD(j*2)
        !         open(14, file='def_coords2.txt')
        !         write(14,'(F36.16, ",", F36.16)') YD(2*j-1), YD(2*j)
        !     enddo
        ! enddo
        ! do i = 1,2,1
        !   YD(2*i-1) = 120d0 + RNDM()*(X_size-240d0)
        !   YD(2*i) = Y_size*RNDM()
        !   write(*, '("X", I3, " = ", F7.1, ", Y", I3, " = ", F7.1)') i, YD(i*2-1), i, YD(i*2)
        !   write(12,'(F36.16, ",", F36.16)') YD(2*i-1), YD(2*i)
        ! enddo
        !Triangular distribution
        count = 0
        do i1 = 1, 36, 1
          yy = 120d0*(i1-1) + 60d0 ! снизу вверх заполнение 
          count = count + 1
          if (mod(i1+10,2) /=0 ) then
            m = 28
          else
            m = 27
          endif
          t1 = ((count-1)/2)*27 + ((count-1) - (count-1)/2)*28 + 1
          !print*, t1
          count1 = 0
          do j1 = t1, t1 + m-1, 1
              count1 = count1 + 1
              if (mod(i1+10,2) /= 0) then
                YD(j1*2-1) = 120d0*(count1-1) + 540d0     
              else
                YD(j1*2-1) = 120d0*(count1-1) + 600d0
              endif
              YD(2*j1) = yy
              write(*, '("X", I3, " = ", F7.1, ", Y", I3, " = ", F7.1)') j1, YD(j1*2-1), j1, YD(j1*2)
              open(13, file='def_coords1.txt')
              write(13,'(F36.16, ",", F36.16)') YD(2*j1-1), YD(2*j1)
          enddo
      enddo
      return
    end subroutine init_config

    subroutine add_vortices(system_energy)
      implicit none
      integer::delta_N, M_old1, M_new, i, j, i2
      real*8::delta_H, x_paart, y_paart, diss, system_energy
      real*8, allocatable:: Y_newvort(:)
      system_energy = 0d0
      delta_H = 100d0
      delta_N = 90
      M_old1 = size(Y)
      allocate(Y_newvort(M_old1))
      Y_newvort = Y
      deallocate(Y)
      Pcle_Num = Pcle_Num + delta_N
      allocate(Y(2*Pcle_Num))
      M_new = 2*Pcle_Num
      do i2 = 1, M_old1, 1
        Y(i2) = Y_newvort(i2)
      enddo
      do i = M_old1+1, M_new-1, 2        
        if (i <= delta_N + M_old1 + 1) then
          x_paart = 120d0* RNDM()       
        else 
          x_paart = 120d0* (RNDM()-1) + X_size
        endif
        y_paart = Y_size * RNDM()
        Y(i) = x_paart         
        Y(i+1) = y_paart  
        do j = i-2, 1, -2
            call distance(x_paart, y_paart, Y(j), Y(j+1), diss)
            if ((x_paart == Y(j)  .AND.  y_paart == Y(j+1)) .OR.&
                (diss == 2)) then       
                print*, "Rechoosing coordinates"              
                x_paart = 120d0 * RNDM()
                Y(i) = x_paart                 
            endif
        end do
        write(*, '("X", I4, " = ", F7.1, ", Y", I4, " = ", F7.1)') i, x_paart, i, y_paart
      end do
    end subroutine add_vortices

    subroutine delete_vortices(system_energy)
      implicit none
      integer :: i, M, j, count_in, i1
      real*8 :: system_energy
      real*8, allocatable :: Y_(:)
      system_energy = 0
      M = size(Y)
      count_in = 0
      do i = 1, M-1, 2
        if (Y(i) <= 120d0 .OR. Y(i) >= (X_size - 120d0)) then
          do j = i, M-1, 1
            Y(i) = Y(i+1)
          end do
          count_in = count_in + 1
        end if 
      end do
      allocate(Y_(M-count_in))
      do i1 = 1, M-count_in, 1
        Y_(i1) = Y(i1)
      end do
      deallocate(Y)
      allocate(Y(M-count_in))
      Y = Y_
    end subroutine delete_vortices

    subroutine B_calc(Ya, suma)
      implicit none
      real*8, intent(inout) :: Ya(:)
      real*8 :: x_chst, y_chst, suma, drk, k_chst
      integer :: M_all, i, j, k, Mx, My, N
      M_all = size(Ya)
      Mx = X_size/2
      My = Y_size/2
      N = Mx*My
      suma = 0d0
      do i = 1, Mx, 1
        x_chst = 1d0 + (i-1)*2d0
        do j = 1, My, 1
          y_chst = 1d0 + (j-1)*2d0
          do k = 1, M_all-1, 2
            call distance(Y(k), Y(k+1), x_chst, y_chst, drk)
            call ik01a(lambda_inv*drk, k_chst)
            suma = suma + k_chst
          end do
        end do
      end do
      suma = suma/N
    end subroutine B_calc


    subroutine energy_calc(NumOf_pcle1, system_energy1) 
        implicit none                                      
        integer :: NumOf_pcle1, j2, j3                        
        real*8 :: system_energy1, dr1
        system_energy1 = 0d0 
        !For Lennard-Jons
        do j2 = 1, 4*NumOf_pcle1-6, 4
            do j3 = j2 + 2, 4*NumOf_pcle1-3
                call distance(Y(j2), Y(j2+1), &
                              Y(j3), Y(j3+1), dr1)
                system_energy1 = system_energy1 + 4*((sigma**12)/(dr1**12) - (sigma**6)/(dr1**6))
            end do
        end do
        !For vortices
        do j2 = 1, 2*NumOf_pcle1-3, 2
          do j3 = j2 + 2, 2*NumOf_pcle1-1
              call distance(Y(j2), Y(j2+1), &
                            Y(j3), Y(j3+1), dr1)
              system_energy1 = system_energy1 + 4*((sigma**12)/(dr1**12) - (sigma**6)/(dr1**6))
          end do
        end do
        return
    end subroutine energy_calc 
    
    subroutine distance(x_1, y_1, x_2, y_2, r)
        implicit none
        real*8 :: x_1, x_2, y_1, y_2
        real*8 :: r
        !print*, x_1, y_2
        r = dsqrt(((x_2-x_1)**2)+((y_2-y_1)**2))
        return
      end subroutine distance

      subroutine interaction_LJ(ind, Y3,  DY_xpart, DY_ypart)
        !считает взаимодействие одной выбранной частицы с индексом ind 
        !со всеми остальными
        implicit none
        integer :: ind, j3, M3
        real*8 :: DY_xpart, DY_ypart, dr1
        real*8,intent(in) :: Y3(:)
        DY_xpart = 0d0
        DY_ypart = 0d0
        M3 = size(Y3)
        do j3 = 1, M3-3, 4
          if (j3 == ind) cycle
          call distance(Y3(j3), Y3(j3+1), Y3(ind), Y3(ind+1), dr1)
          if(dr1 < 0.1d0) dr1 = 0.1d0 !попробуй без обрезки, было 0.35
            DY_xpart = DY_xpart + 8  * esh * (Y3(ind) - Y3(j3)) * (6 * (sigma**12) /(dr1**14) - 3 * (sigma**6)/(dr1**8)) !+
            DY_ypart = DY_ypart + 8 * esh * (Y3(ind+1) - Y3(j3+1)) * (6 *(sigma**12)/(dr1**14) - 3 * (sigma**6)/(dr1**8))
        end do
        return
      end subroutine interaction_LJ

      subroutine interaction_vortices(ind, Y3,  DY_xpart, DY_ypart)
        !считает взаимодействие одной выбранной частицы с индексом ind 
        !со всеми остальными
        implicit none
        integer :: ind, j3, M3
        real*8 :: DY_xpart, DY_ypart, dr1, k1, dr2, k2, dr3, k3, dr4, dr5, k4, k5
        real*8,intent(inout) :: Y3(:)
        real*8 :: dr6, k6, dr7, k7, dr8, k8, dr9, k9
        DY_xpart = 0d0
        DY_ypart = 0d0
        M3 = size(Y3)
        do j3 = 1, M3-1, 2
          call distance(Y3(j3), Y3(j3+1), Y3(ind), Y3(ind+1), dr1)
          if (dr1 < 2.2d0) then
            dr1 = 2.2d0
          endif
          if (dr1 <= 720d0 .AND. j3 /= ind) then
            call ik01a(lambda_inv*dr1, k1)
            DY_xpart = DY_xpart + (k1 * (Y3(ind) - Y3(j3)))/dr1 
            DY_ypart = DY_ypart + (k1 * (Y3(ind+1) - Y3(j3+1)))/dr1
          end if 
          call distance(Y3(j3), Y3(j3+1)-Y_size, Y3(ind), Y3(ind+1), dr2)
          if (dr2 <= Y_size*0.5d0) then
            call ik01a(lambda_inv*dr2, k2)
            DY_xpart = DY_xpart + (k2 * (Y3(ind) - Y3(j3)))/dr2 
            DY_ypart = DY_ypart + (k2 * (Y3(ind+1) - (Y3(j3+1)-Y_size)))/dr2
          end if
          call distance(Y3(j3), Y3(j3+1)+Y_size, Y3(ind), Y3(ind+1), dr3)
          if (dr3 <= Y_size*0.5d0) then
            call ik01a(lambda_inv*dr3, k3)
            DY_xpart = DY_xpart + (k3 * (Y3(ind) - Y3(j3)))/dr3 
            DY_ypart = DY_ypart + (k3 * (Y3(ind+1) - (Y3(j3+1)+Y_size)))/dr3
          end if
          call distance(Y3(j3)+X_size, Y3(j3+1), Y3(ind), Y3(ind+1), dr4)
          if (dr4 <= X_size*0.5d0) then
            call ik01a(lambda_inv*dr4, k4)
            DY_xpart = DY_xpart + (k4 * (Y3(ind) - (Y3(j3)+X_size)))/dr4 
            DY_ypart = DY_ypart + (k4 * (Y3(ind+1) - Y3(j3+1)))/dr4
          end if
          call distance(Y3(j3)-X_size, Y3(j3+1), Y3(ind), Y3(ind+1), dr5)
          if (dr5 <= X_size*0.5d0) then
            call ik01a(lambda_inv*dr5, k5)
            DY_xpart = DY_xpart + (k5 * (Y3(ind) - (Y3(j3)-X_size)))/dr5
            DY_ypart = DY_ypart + (k5 * (Y3(ind+1) - Y3(j3+1)))/dr5
          end if
          call distance(Y3(j3)-X_size, Y3(j3+1)-Y_size, Y3(ind), Y3(ind+1), dr6)
          if (dr6 <= X_size*0.5d0) then
            call ik01a(lambda_inv*dr6, k6)
            DY_xpart = DY_xpart + (k6 * (Y3(ind) - (Y3(j3)-X_size)))/dr6 
            DY_ypart = DY_ypart + (k6 * (Y3(ind+1) - (Y3(j3+1)-Y_size)))/dr6
          end if
          call distance(Y3(j3)+X_size, Y3(j3+1)+Y_size, Y3(ind), Y3(ind+1), dr7)
          if (dr7 <= X_size*0.5d0) then
            call ik01a(lambda_inv*dr7, k7)
            DY_xpart = DY_xpart + (k7 * (Y3(ind) - (Y3(j3)+X_size)))/dr7
            DY_ypart = DY_ypart + (k7 * (Y3(ind+1) - (Y3(j3+1)+Y_size)))/dr7
          end if
          call distance(Y3(j3)-X_size, Y3(j3+1)+Y_size, Y3(ind), Y3(ind+1), dr8)
          if (dr8 <= X_size*0.5d0) then
            call ik01a(lambda_inv*dr8, k8)
            DY_xpart = DY_xpart + (k8 * (Y3(ind) - (Y3(j3)-X_size)))/dr8 
            DY_ypart = DY_ypart + (k8 * (Y3(ind+1) - (Y3(j3+1)+Y_size)))/dr8
          end if
          call distance(Y3(j3)+X_size, Y3(j3+1)-Y_size, Y3(ind), Y3(ind+1), dr9)
          if (dr9 <= X_size*0.5d0) then
            call ik01a(lambda_inv*dr9, k9)
            DY_xpart = DY_xpart + (k9 * (Y3(ind) - (Y3(j3)+X_size)))/dr9
            DY_ypart = DY_ypart + (k9 * (Y3(ind+1) - (Y3(j3+1)-Y_size)))/dr9
          end if
          !Просто расчитываем силу взаимодействия вихря с ближайшими соседями, которые как бы в верхней ячейке
          !и в нижней ячейке. Распределение вихрей вверху и внизу изменяется аналогично изменениям на основной ячейке
          !(основная - та которую видим на графиках). Поэтому туда и дефекты не нужно забрасывать!!!
        end do
        return
      end subroutine interaction_vortices

      subroutine interaction_surf_vort(ind, Y3, DY_xpart, DY_ypart)
        implicit none
        integer :: ind, j3, M3
        real*8 :: DY_xpart, DY_ypart, dr1, k1, k2, k_x, k_d_x, dr2
        real*8,intent(in) :: Y3(:)
        DY_xpart = 0d0
        DY_ypart = 0d0
        M3 = size(Y3)
        call ik01a(2d0*lambda_inv*Y3(ind), k_x) !взаимодействие с тем, кто за нулевой границей
        call ik01a(2d0*lambda_inv*(X_size-Y3(ind)), k_d_x) !взаимодействие с тем, кто за X_size
        DY_xpart = k_x + k_d_x
        do j3 = 1, M3-1, 2
          if (j3 == ind) cycle
          call distance(2*X_size-Y3(j3), Y3(j3+1), Y3(ind), Y3(ind+1), dr1) !взимодействие с тем, кто за X_size
          call distance(0d0, Y3(j3+1), Y3(ind)+Y3(j3), Y3(ind+1), dr2) !взимодействие с тем, кто за нулевой границей
          if (dr1 <= 2.2d0) dr1 = 2.2d0
          if (dr2 <= 2.2d0) dr2 = 2.2d0
          if (dr1 <= 720d0) then           !Перепроверь здесь, отсюда Nan идет
            call ik01a(lambda_inv*dr1, k1)
          else
            k1 = 0d0
          endif
          if (dr2 <= 720d0) then
            call ik01a(lambda_inv*dr2, k2)
          else
            k2 = 0d0
          endif
          DY_xpart = DY_xpart + (k1 * (Y(ind)-2*X_size+Y3(j3)))/dr1  
          DY_ypart = DY_ypart + (k1 * (Y(ind+1) - Y(j3+1)))/dr1
        end do
        return
     end subroutine interaction_surf_vort

      subroutine int_wth_defects(ind1, YR, X_part, Y_part)
        implicit none
        integer :: ind1, Numb_def,i 
        real*8 :: X_part, Y_part
        real*8 :: dr1, ksi_p, ksi_p_inv
        real*8,intent(inout) :: YR(:)
        ksi_p = 25d0 !нм
        ksi_p_inv = 1/ksi_p
        X_part = 0d0
        Y_part = 0d0
        ! просто координаты дефектов
        Numb_def = 990
        do i = 1, 2*Numb_def-1, 2
            call distance(YR(ind1), YR(ind1+1), YD(i), YD(i+1), dr1) 
            if (dr1 <= ksi_p) then
                X_part = X_part + (YR(ind1) - YD(i))*ksi_p_inv 
                Y_part = X_part + (YR(ind1+1) - YD(i+1))*ksi_p_inv 
            endif
        enddo
        !здесь сразу считается взаимодействие всех частиц с Numb_def дефектами
        !y часть градиента и х часть градиента
        return
    end subroutine int_wth_defects




    subroutine ik01a(x,  bk1)
        implicit none

        real*8, save, dimension (12) :: a = (/ &
        0.125D+00, 7.03125D-02, &
        7.32421875D-02, 1.1215209960938D-01, &
        2.2710800170898D-01, 5.7250142097473D-01, &
        1.7277275025845D+00, 6.0740420012735D+00, &
        2.4380529699556D+01, 1.1001714026925D+02, &
        5.5133589612202D+02, 3.0380905109224D+03 /)
        real*8, save, dimension (8) :: a1 = (/ &
        0.125D+00, 0.2109375D+00, &
        1.0986328125D+00, 1.1775970458984D+01, &
        2.1461706161499D+02, 5.9511522710323D+03, &
        2.3347645606175D+05, 1.2312234987631D+07 /)
        real*8, save, dimension (12) :: b = (/ &
        -0.375D+00, -1.171875D-01, &
        -1.025390625D-01, -1.4419555664063D-01, &
        -2.7757644653320D-01, -6.7659258842468D-01, &
        -1.9935317337513D+00, -6.8839142681099D+00, &
        -2.7248827311269D+01, -1.2159789187654D+02, &
        -6.0384407670507D+02, -3.3022722944809D+03 /)
        real*8 :: bi0, bi1, bk0, bk1, ca, cb, ct, di0, di1, dk0, dk1, el
        integer*4 :: k, k00
        real*8 :: P, r, w0, ww, x, x2, xr, xr2

        P = 3.1415926535897932D+00
        el = 0.5772156649015329D+00
        x2 = x * x
        if (x == 0.0D+00) then
            bi0 = 1.0D+00
            bi1 = 0.0D+00
            bk0 = 1.0D+300
            bk1 = 1.0D+300
            di0 = 0.0D+00
            di1 = 0.5D+00
            dk0 = -1.0D+300
            dk1 = -1.0D+300
            return
        elseif (x <= 18.0D+00) then
            bi0 = 1.0D+00
            r = 1.0D+00
            do k = 1, 50
                r = 0.25D+00 * r * x2 / (k * k)
                bi0 = bi0 + r
                if (abs (r / bi0) < 1.0D-15) exit
            end do
            bi1 = 1.0D+00
            r = 1.0D+00
            do k = 1, 50
                r = 0.25D+00 * r * x2 / (k * (k + 1))
                bi1 = bi1 + r
                if (abs (r / bi1) < 1.0D-15) exit
            end do
            bi1 = 0.5D+00 * x * bi1
        else
            if (x < 35.0D+00) then
                k00 = 12
            elseif (x < 50.0D+00) then
                k00 = 9
            else
                k00 = 7
            end if
            ca = dexp (x) / dsqrt (2.0D+00 * P * x)
            bi0 = 1.0D+00
            xr = 1.0D+00 / x
            do k = 1, k00
                bi0 = bi0 + a(k) * xr ** k
            end do
            bi0 = ca * bi0
            bi1 = 1.0D+00
            do k = 1, k00
                bi1 = bi1 + b(k) * xr ** k
            end do
            bi1 = ca * bi1
        end if

        if (x <= 9.0D+00) then
            ct = - (dlog (x / 2.0D+00) + el)
            bk0 = 0.0D+00
            w0 = 0.0D+00
            r = 1.0D+00
            do k = 1, 50
                w0 = w0 + 1.0D+00 / k
                r = 0.25D+00 * r / (k * k) * x2
                bk0 = bk0 + r * (w0 + ct)
                if (abs ((bk0 - ww) / bk0) < 1.0D-15) exit
                ww = bk0
            end do
            bk0 = bk0 + ct
        else
            cb = 0.5D+00 / x
            xr2 = 1.0D+00 / x2
            bk0 = 1.0D+00
            do k = 1, 8
                bk0 = bk0 + a1(k) * xr2 ** k
            end do
            bk0 = cb * bk0 / bi0
        end if

        bk1 = (1.0D+00 / x - bi1 * bk0) / bi0
        di0 = bi1
        di1 = bi0 - bi1 / x
        dk0 = - bk1
        dk1 = - bk0 - bk1 / x

        return
    end subroutine ik01a

     
end module processes1