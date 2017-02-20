      subroutine FCN(T,Z,DZ,TRPAR)
      doubleprecision T
      doubleprecision Z(3)
      doubleprecision DZ(3)
      doubleprecision TRPAR(2)

        DZ(1) = TRPAR(1)*exp(0.1D0*Z(3))-Z(1)*Z(2)**2-0.55D-2*Z(1)
        DZ(2) = Z(1)*Z(2)**2+0.55D-2*Z(1)-Z(2)
        DZ(3) = Z(2)-TRPAR(2)*Z(3)
        return
      end
