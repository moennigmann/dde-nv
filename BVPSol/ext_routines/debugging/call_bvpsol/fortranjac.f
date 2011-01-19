      subroutine DFDY (T, Z, DF, TRPAR)
        doubleprecision T
        doubleprecision Z(4)
        doubleprecision DF(4,4)
        doubleprecision TRPAR(2)
        DF(1, 1) = -0.1D0 - TRPAR(2) * Z(2) * Z(4)
        DF(1, 2) = -TRPAR(2) * Z(1) * Z(4)
        DF(1, 3) = 0
        DF(1, 4) = -TRPAR(2) * Z(1) * Z(2)
        DF(2, 1) = -TRPAR(2) * Z(2) * Z(4)
        DF(2, 2) = -TRPAR(1) * Z(3) - TRPAR(2) * Z(1) * Z(4)
        DF(2, 3) = -TRPAR(1) * Z(2)
        DF(2, 4) = -TRPAR(2) * Z(1) * Z(2)
        DF(3, 1) = 0.3D1 * TRPAR(2) * Z(2) * Z(4)
        DF(3, 2) = TRPAR(1) * Z(3) + 0.3D1 * TRPAR(2) * Z(1) * Z(4)
        DF(3, 3) = TRPAR(1) * Z(2) - 0.1000D4 * Z(3) - 0.20D2
        DF(3, 4) = 0.3D1 * TRPAR(2) * Z(1) * Z(2)
        DF(4, 1) = -TRPAR(2) * Z(2) * Z(4)
        DF(4, 2) = -TRPAR(2) * Z(1) * Z(4)
        DF(4, 3) = 0.1000D4 * Z(3)
        DF(4, 4) = -0.535D1 - TRPAR(2) * Z(1) * Z(2)
        return
      end
