
This folder contains some useful tools for developers:

  create_methods.php -  Creates methods for a class:
                        member accessors, default constructors, destructor, assignment operator and equality operators.

  check_test         -  Base tool for checking the tests. It can be built by executing 'make'.

  create_test.php    -  Creates a test for a file/class using 'check_test'.

  correct_test.php   -  Helps the user to correct CHECK macros of tests using 'check_test'.
	                      Function names inside the CHECK macros are automatically replaced by
												the correct declarations in the header file.
  
  check_includes.php -  Checks for unneeded includes.

  checker.php        -  Reports errors in the code, test, documentation, ...

Unfortunately there are some bugs in the C++ parser 'check_test'.
It fails under the following conditions:

  - Extra semicolon after method implementation in the header.
  
  - Method arguments are distributed over several rows.
  
  - Nested class declarations.
	
	- Several variable declarations in one row (int a,b,c;)
	
	- Union declarations

	- Pure access declarations

