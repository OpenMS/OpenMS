Adding pure python classes/functionality
========================================

Pure python modules can be found `here <https://github.com/OpenMS/OpenMS/tree/develop/src/pyOpenMS/pyopenms>`.
Just add a new `.py` file and a new submodule will appear that can be imported with `import pyopenms.submodule`.

Remember to add additional requirements in the setup.py but try to avoid large dependencies
unless absolutely necessary.

Testing is done via pytest.
Every new module, class, function, member should be documented with Sphinx reStructuredText
docstrings. See the `Sphinx-RTD-Tutorial <https://sphinx-rtd-tutorial.readthedocs.io/en/latest/docstrings.html>` and the `Python Developers Guide <https://devguide.python.org/documentation/start-documenting/index.html>`.
