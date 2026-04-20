Help us to make OpenMS better and become part of the OpenMS open-source community.

This document is displayed because you either opened an issue or you want to provide your code as a pull request for inclusion into OpenMS. Please take a look at the appropriate section below to find some details on how we handle this process.
When interacting with other developers, users or anyone else from our community, please adhere to
[our CODE OF CONDUCT](https://github.com/OpenMS/OpenMS/blob/develop/CODE_OF_CONDUCT.md)

> **New to open source or the project?** Start with our [Contributor Onboarding Guide](https://openms.de/contribute/) for a beginner-friendly introduction.

# Reporting an Issue:

You most likely came here to:
  - report bugs or annoyances
  - pose questions
  - point out missing documentation
  - request new features

To do so, open a [new issue](https://github.com/OpenMS/OpenMS/issues/new/choose) describing the bug, your question, etc.

If you found a bug, e.g. an OpenMS tool crashes during data processing, it is essential to provide some basic information:
  - the OpenMS version you are running
  - the platform you are running OpenMS on (Windows 10, ...)
  - how you installed OpenMS (e.g., from within KNIME, binary installers, self compiled)
  - a description on how to reproduce the bug
  - relevant tool output (e.g., error messages)
  - data to repoduce the bug (If possible as a GitHub gist. Other platforms like Dropbox, Google Drive links also work. If you can't share the data publicly please indicate this and we will contact you in private.)

If you are an official OpenMS team member:
  - label your issue using github labels (e.g. as: question, defect) that indicate the type of issue and which components of OpenMS (blue labels) are affected. The severity is usually assigned by OpenMS maintainers and used internally to e.g. indicate if a bug is a blocker for a new release.

# Opening a Pull Request

Before getting started we recommend taking a look at the developer guide:
https://openms.readthedocs.io/en/latest/manual/develop.html

Before you open the pull request, make sure you
 - adhere to [our coding conventions](https://openms.readthedocs.io/en/latest/manual/develop.html#coding-conventions)
 - have [unit tests and functional tests](https://openms.readthedocs.io/en/latest/manual/develop.html#automated-unit-tests) - see also [this example](https://github.com/OpenMS/OpenMS/blob/develop/src/tests/class_tests/openms/source/MSNumpressCoder_test.cpp)
 - have [proper documentation](https://openms.readthedocs.io/en/latest/manual/develop/developer-faq.html#doxygen-documentation) - see also [this example](https://github.com/OpenMS/OpenMS/blob/develop/src/openms/include/OpenMS/FORMAT/MSNumpressCoder.h)
 - have Python bindings — nanobind binding files are in `src/pyOpenMS/bindings/`; see `src/pyOpenMS/CLAUDE.md` for wrapping instructions
A core developer will review your changes to the main development branch (develop) and approve them (or ask for modifications). You may indicate the prefered reviewer(s) by adding links to them in a comment section (e.g., @cbielow @hendrikweisser @hroest @jpfeuffer @timosachsenberg)

Also consider getting in contact with the core developers early. They might provide additional guidance and valuable information on how your specific aim is achieved. This might give you a head start in, for example, developing novel tools or algorithms.

Happy coding!

---

## Additional Resources

For detailed coding conventions, architectural guidelines, and technical documentation, see the [OpenMS Developer Guide](https://openms.readthedocs.io/en/latest/contribute-to-openms/index.html).
