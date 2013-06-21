      subroutine DFDY(T,Z,DF,TRPAR)
      doubleprecision T
      doubleprecision Z(3)
      doubleprecision DF(3,3)
      doubleprecision TRPAR(2)

        DF(1,1) = -Z(2)**2-0.55D-2
        DF(1,2) = -2*Z(1)*Z(2)
        DF(1,3) = 0.1D0*TRPAR(1)*exp(0.1D0*Z(3))
        DF(2,1) = Z(2)**2+0.55D-2
        DF(2,2) = 2*Z(1)*Z(2)-1
        DF(2,3) = 0
        DF(3,1) = 0
        DF(3,2) = 1
        DF(3,3) = -TRPAR(2)
        return
      end
