      subroutine funcJ0P(Z,DJP,TRPAR)
      doubleprecision Z(3)
      doubleprecision DJP(3,3,2)
      doubleprecision TRPAR(2)

        DJP(1,1,1) = 0
        DJP(1,1,2) = 0
        DJP(1,2,1) = 0
        DJP(1,2,2) = 0
        DJP(1,3,1) = 0.1D0*exp(0.1D0*Z(3))
        DJP(1,3,2) = 0
        DJP(2,1,1) = 0
        DJP(2,1,2) = 0
        DJP(2,2,1) = 0
        DJP(2,2,2) = 0
        DJP(2,3,1) = 0
        DJP(2,3,2) = 0
        DJP(3,1,1) = 0
        DJP(3,1,2) = 0
        DJP(3,2,1) = 0
        DJP(3,2,2) = 0
        DJP(3,3,1) = 0
        DJP(3,3,2) = -1
        return
      end
