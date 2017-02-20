      subroutine FCN (T, Z, DZ, TRPAR)
        doubleprecision T
        doubleprecision Z(4)
        doubleprecision DZ(4)
        doubleprecision TRPAR(2)
        DZ(1) = 0.8D0 - 0.1D0 * Z(1) - TRPAR(2) * Z(1) * Z(2) * Z(4)
        DZ(2) = 0.825D0 - TRPAR(1) * Z(2) * Z(3) - TRPAR(2) * Z(1) * Z(2
     #) * Z(4)
        DZ(3) = TRPAR(1) * Z(2) * Z(3) - 0.500D3 * Z(3) ** 2 + 0.3D1 * T
     #RPAR(2) * Z(1) * Z(2) * Z(4) - 0.20D2 * Z(3) + 0.1D1 / 0.100000D6
        DZ(4) = 0.500D3 * Z(3) ** 2 - 0.535D1 * Z(4) - TRPAR(2) * Z(1) *
     # Z(2) * Z(4)
        return
      end
