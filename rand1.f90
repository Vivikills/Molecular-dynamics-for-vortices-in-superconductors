module rand
    implicit none
    
    integer :: i97, j97
    double precision, parameter :: r1 = 5.d-15, r2 = 1.d-14
    real*8 :: c, cd, cm, Ugen(97), rval
    
    
    contains
     !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     ! RANDOMIZER SEED ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        subroutine SRANMAR(IJ,  KL)
            implicit none
            
            integer :: ij, kl, i7, k7, j7, l7, ii7, jj7, m7
            real*8  :: s, tt
            
            if (ij < 0   .OR.   ij > 31328   .OR.   &
                    kl < 0   .OR.   kl > 30081) then
                print '(A)', " The first  seed must be between 0 and 31328"
                print '(A)', " The second seed must be between 0 and 30081"
                STOP
            end if
            
            i7 = MOD(ij/177, 177) + 2
            j7 = MOD(ij, 177) + 2
            k7 = MOD(kl/169, 178) + 1
            l7 = MOD(kl, 169)
            
            do ii7 = 1, 97
                s = 0d0
                tt = 0.5d0
                do jj7 = 1, 24
                    m7 = MOD(mod(i7*j7, 179)*k7, 179)
                    i7 = j7
                    j7 = k7
                    k7 = m7
                    l7 = MOD(53*l7+1, 169)
                    if (MOD(l7*m7, 64) > 32) s = s + tt
                    tt = 0.5d0 * tt
                end do
                UGEN(ii7) = s
            end do
            
            c  = 362436d0   / 16777216d0
            cd = 7654321d0  / 16777216d0
            cm = 16777213d0 / 16777216d0
            i97 = 97
            j97 = 33
            
            RETURN
        end subroutine SRANMAR
        
     !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     ! RANDOMIZER FUNCTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        double precision function RNDM()
            implicit none
            
            rval = ugen(i97) - ugen(j97)
            if (rval < 0d0) rval = rval + 1d0
            ugen(i97) = rval
            i97 = i97 - 1
            if (i97 == 0) i97 = 97            
            j97 = j97 - 1
            if (j97 == 0) j97 = 97
            c = c - cd            
            if (c < 0d0) c = c + cm
            rval = rval - c            
            if (rval < 0d0) rval = rval + 1d0
            RNDM = max(rval-r1, r2)
            
            RETURN
        end function RNDM
end module rand

