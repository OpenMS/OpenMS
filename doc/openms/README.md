OpenMS Documentation
====================

This repository is the home for OpenMS documentation. OpenMS documentation is categorised as per:

- Supported libraries, tools, packages.
- Users who prefer GUI.
- Users who prefer working with code.
- Examples and tutorials of working with OpenMS.
- Adding functionalities to OpenMS.
- Extensibility of OpenMS.
- Developer and contributor documentation.
- Workshop/Training tutorial.
- Quick start user and developer guides.
- Most efficient and user centric OpenMS tools.
- and a lot of other interesting reads with sample data!

The documentation can be browsed online and offline.

## Browse OpenMS documentation online

1. https://openms.readthedocs.io/en/develop
2. https://openms.readthedocs.io/en/latest

`develop` denotes the developing version of OpenMS documentation. `latest` corresponds to the OpenMS stable release
documentation.


## Browse OpenMS documentation offline

1. Download latest [PDF](https://openms.readthedocs.io/_/downloads/en/latest/pdf/).
2. Download latest [HTMLZip](https://openms.readthedocs.io/_/downloads/en/latest/htmlzip/).
3. Download develop [PDF](https://openms.readthedocs.io/_/downloads/en/develop/pdf/).
4. Download develop [HTMLZip](https://openms.readthedocs.io/_/downloads/en/develop/htmlzip/).

## Build OpenMS Docs locally

In the root directory, run

```
pip3 install -r requirements.txt
```

to install dependencies. It is suggested to create a virtual environment, for more information read [venv](https://docs.python.org/3/library/venv.html) docs.

As the next step, there are a few options to build OpenMS/OpenMS-docs locally:

1. Run `make html` in root directory.
2. Run `sphinx-build -a docs /tmp/rtd` to build all files.
3. Run `sphinx-build docs/ /tmp/rtd` to build changed files.

Now, open index.html file to view the changes or output of OpenMS ReadTheDocs documentation.

## Contributing to OpenMS Docs

Please read our [contributing guidelines](.github/CONTRIBUTING.md), before starting with contributing to OpenMS
Documentation.

Let us know what you would like to read in OpenMS documentation using [GitHub issues](https://github.com/OpenMS/OpenMS-docs/issues/new/choose)!
