      subroutine DGSUB(ICount,NCOMP,Z,DG,TRPAR,IPAR)
      integer ICount
      integer NCOMP
      doubleprecision Z(6)
      doubleprecision DG(6)
      doubleprecision TRPAR(2)
      integer IPAR

      doubleprecision DGSET(6,6)
      integer ind

        DGSET(1,1) = 0
        DGSET(1,2) = 1
        DGSET(1,3) = 0
        DGSET(1,4) = 0
        DGSET(1,5) = 0
        DGSET(1,6) = 0
        DGSET(2,1) = 1
        DGSET(2,2) = 0
        DGSET(2,3) = 0
        DGSET(2,4) = 0
        DGSET(2,5) = -1
        DGSET(2,6) = 0
        DGSET(3,1) = 0
        DGSET(3,2) = 0
        DGSET(3,3) = 1
        DGSET(3,4) = 0
        DGSET(3,5) = 0
        DGSET(3,6) = -1
        DGSET(4,1) = 0
        DGSET(4,2) = 1
        DGSET(4,3) = 0
        DGSET(4,4) = 0
        DGSET(4,5) = 0
        DGSET(4,6) = 0
        DGSET(5,1) = 1
        DGSET(5,2) = 0
        DGSET(5,3) = 0
        DGSET(5,4) = 0
        DGSET(5,5) = -1
        DGSET(5,6) = 0
        DGSET(6,1) = 0
        DGSET(6,2) = 0
        DGSET(6,3) = 1
        DGSET(6,4) = 0
        DGSET(6,5) = 0
        DGSET(6,6) = -1
        do 1000 ind = 1,6,1
          DG(ind) = DGSET(ICount,ind)
1000    continue
        return
      end
