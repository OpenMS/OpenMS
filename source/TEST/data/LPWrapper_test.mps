NAME          LPWrapper_test.lp
ROWS
 N  OBJROW  
 L  Const_1 
 L  Const_2 
 L  Const_3 
COLUMNS
    MARK0000  'MARKER'                 'INTORG'
    y         OBJROW                          1
    y         Const_1                         1
    y         Const_2                         2
    y         Const_3                         3
    x         Const_1                        -1
    x         Const_2                         3
    x         Const_3                         2
    MARK0001  'MARKER'                 'INTEND'
RHS
    rhs       Const_2                        12
    rhs       Const_3                        12
BOUNDS
 LI bnd       y         0
 LI bnd       x         0
ENDATA
