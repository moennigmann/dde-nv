      subroutine DFSUB(NCOMP,X,Z,DF,TRPAR,IPAR)
      integer NCOMP
      doubleprecision X
      doubleprecision Z(6)
      doubleprecision DF(6,6)
      doubleprecision TRPAR(2)
      integer IPAR

        DF(1,1) = Z(4)*(-Z(2)**2-0.55D-2)
        DF(1,2) = -2*Z(4)*Z(1)*Z(2)
        DF(1,3) = 0.1D0*Z(4)*TRPAR(1)*exp(0.1D0*Z(3))
        DF(1,4) = TRPAR(1)*exp(0.1D0*Z(3))-Z(1)*Z(2)**2-0.55D-2*Z(1)
        DF(1,5) = 0
        DF(1,6) = 0
        DF(2,1) = Z(4)*(Z(2)**2+0.55D-2)
        DF(2,2) = Z(4)*(2*Z(1)*Z(2)-1)
        DF(2,3) = 0
        DF(2,4) = Z(1)*Z(2)**2+0.55D-2*Z(1)-Z(2)
        DF(2,5) = 0
        DF(2,6) = 0
        DF(3,1) = 0
        DF(3,2) = Z(4)
        DF(3,3) = -Z(4)*TRPAR(2)
        DF(3,4) = Z(2)-TRPAR(2)*Z(3)
        DF(3,5) = 0
        DF(3,6) = 0
        DF(4,1) = 0
        DF(4,2) = 0
        DF(4,3) = 0
        DF(4,4) = 0
        DF(4,5) = 0
        DF(4,6) = 0
        DF(5,1) = 0
        DF(5,2) = 0
        DF(5,3) = 0
        DF(5,4) = 0
        DF(5,5) = 0
        DF(5,6) = 0
        DF(6,1) = 0
        DF(6,2) = 0
        DF(6,3) = 0
        DF(6,4) = 0
        DF(6,5) = 0
        DF(6,6) = 0
        return
      end
