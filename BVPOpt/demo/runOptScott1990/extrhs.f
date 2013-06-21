      subroutine FSUB(NCOMP,X,Z,F,TRPAR,IPAR)
      integer NCOMP
      doubleprecision X
      doubleprecision Z(6)
      doubleprecision F(6)
      doubleprecision TRPAR(2)
      integer IPAR

        F(1) = Z(4)*(TRPAR(1)*exp(0.1D0*Z(3))-Z(1)*Z(2)**2-0.55D-2*Z(1))
        F(2) = Z(4)*(Z(1)*Z(2)**2+0.55D-2*Z(1)-Z(2))
        F(3) = Z(4)*(Z(2)-TRPAR(2)*Z(3))
        F(4) = 0
        F(5) = 0
        F(6) = 0
        return
      end
