
This folder contains some useful tools for developers:

  create_methods.php    -  Creates methods for a class:
                           member accessors, default constructors, destructor, assignment operator and equality operators.

  create_test.php       -  Creates a test for a class.

  correct_test.php      -  Helps the user to correct BEGIN_SECTION macros.
                           Function names inside the BEGIN_SECTION macros are automatically replaced by
                           the correct declarations in the header file.
  
  check_includes.php    -  Checks for unneeded includes.

  checker.php           -  Reports errors in the code, test, documentation, ...
                           This tool relies on doxygen XML output.

  sql2grpah.php         -  Transforms the OpenMS MySQL DB to a GraphML representation.
                           The layout can be done with Yed.
	
  make_dist.sh          -  Is used to create OpenMS release file.

  make_dist_contrib.sh  -  Is used to create the contrib release file.
  
  tools.sh              -  A collection of small bash tools you can 'source' into your shell.


Note: All php scripts are executed by calling them as argument of the php interpreter.
