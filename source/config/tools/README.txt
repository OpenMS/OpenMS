
This folder contains some useful tools for developers:

  create_test.php    -  Creates a test for a file/class using 'check_test'.
  
  missing_tests.php  -  Reports tests that are missing, not executed or probably not complete.
  
  create_methods.php -  Creates methods for a class:
                        member accessors, default constructors, destructor, assignment operator and equality operators.
  
  checker.php        -  Checks if header guards are correct.
                        Checks if tab settings for editors and maintainer is in each file.
                        Checks for unneeded includes.

  correct_test.php   -  Helps the user to correct check headers using 'check_test'.
	
	check_test         -  Base tool for checking the tests. It can be built by executing 'make'.


Unfortunately there are some bugs in the C++ parser 'check_test'.
It fails under the following conditions:

	- Extra semicolon after method implementation in the header.
	
	- Method arguments are distributed over several rows.
	
	- Nested class declarations.