      subroutine GSUB(ICount,NCOMP,Z,G,TRPAR,IPAR)
      integer ICount
      integer NCOMP
      doubleprecision Z(6)
      doubleprecision G
      doubleprecision TRPAR(2)
      integer IPAR

      doubleprecision GSET(6)

        GSET(1) = Z(2)-0.69D-1
        GSET(2) = Z(1)-Z(5)
        GSET(3) = Z(3)-Z(6)
        GSET(4) = Z(2)-0.69D-1
        GSET(5) = Z(1)-Z(5)
        GSET(6) = Z(3)-Z(6)
        G = GSET(ICount)
        return
      end
