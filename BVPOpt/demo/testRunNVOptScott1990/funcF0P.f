      subroutine funcF0P(Z,DFP,TRPAR)
      doubleprecision Z(3)
      doubleprecision DFP(3,2)
      doubleprecision TRPAR(2)

        DFP(1,1) = exp(0.1D0*Z(3))
        DFP(1,2) = 0
        DFP(2,1) = 0
        DFP(2,2) = 0
        DFP(3,1) = 0
        DFP(3,2) = -Z(3)
        return
      end
