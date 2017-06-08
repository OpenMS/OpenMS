NAME          LPWrapper_test.mps
ROWS
 N  OBJROW  
 L  Const_1 
 L  Const_2 
 L  Const_3 
COLUMNS
    y         OBJROW                          1
    y         Const_1                         1
    y         Const_2                         2
    y         Const_3                         3
    x         Const_1                        -1
    x         Const_2                         3
    x         Const_3                         2
RHS
    rhs       Const_2                        12
    rhs       Const_3                        12
BOUNDS
 BV BOUND       y
 BV BOUND       x
ENDATA
