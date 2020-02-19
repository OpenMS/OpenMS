Help us to make OpenMS better and become part of the OpenMS open-source community.

This document is displayed because you either opened an issue or you want to provide your code as a pull request for inclusion into OpenMS. Please take a look at the appropriate section below to find some details on how we handle this process.
When interacting with other developers, users or anyone else from our community, please adhere to
[our CODE OF CONDUCT](https://github.com/OpenMS/OpenMS/blob/develop/CODE_OF_CONDUCT.md)

# Reporting an Issue:

You most likely came here to:
  - report bugs or annoyances
  - pose questions
  - point out missing documentation
  - request new features

If you found a bug, e.g. an OpenMS tool crashes during data processing, it is essential to provide some basic information:
  - the OpenMS version you are running
  - the platform you are running OpenMS on (Windows 10, ...)
  - how you installed OpenMS (e.g., from within KNIME, binary installers, self compiled)
  - a description on how to reproduce the bug
  - relevant tool output (e.g., error messages)
  - data to repoduce the bug (If possible as a GitHub gist. Other platforms like Dropbox, Google Drive links also work. If you can't share the data publicly please indicate this and we will contact you in private.)

If you are an official OpenMS team meber:
  - label your issue using github labels (e.g. as: question, defect) that indicate the type of issue and which components of OpenMS (blue labels) are affected. The severity is usually assigned by OpenMS maintainers and used internally to e.g. indicate if a bug is a blocker for a new release.

# Opening a Pull Request

Before getting started we recommend taking a look at our GitHub-Wiki: https://github.com/OpenMS/OpenMS/wiki#-for-developers

Before you open the pull request, make sure you
 - adhere to [our coding conventions](https://github.com/OpenMS/OpenMS/wiki/Coding-conventions)
 - have [unit tests and functional tests](https://github.com/OpenMS/OpenMS/wiki/Write-tests) - see also [this example](https://github.com/OpenMS/OpenMS/blob/develop/src/tests/class_tests/openms/source/MSNumpressCoder_test.cpp)
 - Have [proper documentation](https://github.com/OpenMS/OpenMS/wiki/Coding-conventions#doxygen) - see also [this examle](https://github.com/OpenMS/OpenMS/blob/develop/src/openms/include/OpenMS/FORMAT/MSNumpressCoder.h) 
 - have [Python bindings](https://pyopenms.readthedocs.io/en/latest/wrap_classes.html) - see also [these instructions](https://github.com/OpenMS/OpenMS/wiki/pyOpenMS#wrap) and [this example](https://github.com/OpenMS/OpenMS/blob/develop/src/pyOpenMS/pxds/MSNumpressCoder.pxd)
 
A core developer will review your changes to the main development branch (develop) and approve them (or ask for modifications). You may indicate the prefered reviewer(s) by adding links to them in a comment section (e.g., @cbielow @hendrikweisser @hroest @jpfeuffer @timosachsenberg)

Also consider getting in contact with the core developers early. They might provide additional guidance and valuable information on how your specific aim is achieved. This might give you a head start in, for example, developing novel tools or algorithms.

Happy coding!
