# Conda

% start-after
## Install via Conda

```{warning}
At this time, we do not provide a conda package for our GUI tools. This means if you want to install e.g., TOPPView or SwathWizard
for use in for example one of our tutorials, please refer to a different installation method below.
```

You can use conda or mamba to install the OpenMS library and tools without user interface. Depending on the conda channel, you can
obtain release versions (`bioconda` channel) and nightly versions (`openms` channel).

1. Follow the instructions to install [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) or [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html).
  In the following, every mention of conda may be substituted by mamba for faster environment solving.

2. We recommend to create a new environment with one of the supported python version versions:
   ```bash
    conda create -n openms python=3.10
   ```

2. Add some channels to find dependencies:
   ```bash
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
   ```
   :::{note}
   You can also add the channels for your current environment only with the `--env` option.
   :::
   :::{warning}
   The order of the channels is important!
   :::
   :::{note}
   conda-forge might already be added if you are using Mambaforge.
   :::


3. Install any of the following packages related to OpenMS

   :::::{tab-set}

   ::::{tab-item} openms
   :sync: openms

   openms contains all OpenMS C++ command-line tools. GUI applications like TOPPView currently cannot be installed via conda.
   ::::

   ::::{tab-item} libopenms
   :sync: libopenms

   libopenms is the C++ library required for the OpenMS C++ Tools to work. This is also an auto-installed dependency of openms.
   ::::

   ::::{tab-item} pyopenms
   :sync: pyopenms

   pyopenms is the python package that allows to use algorithms from libopenms in Python.
   ::::

   ::::{tab-item} openms-thirdparty
   :sync: openms-thirdparty

   openms-thirdparty are external tools that are wrapped in OpenMS with adapters. This package is required
   to use the adapters in the openms package.

   :::{warning}
   Due to unavailability of a large part of the thirdparty tools for macOS via conda, we are not providing
   a openms-thirdparty package on macOS either.
   :::

   ::::

   :::::

   via `bioconda` for release versions


   ::::{tab-set}

   :::{tab-item} openms
   :sync: openms

   ```{code-block} bash 
   conda install openms
   ```
   :::

   :::{tab-item} libopenms
   :sync: libopenms

   ```{code-block} bash
   conda install libopenms
   ```
   :::

   :::{tab-item} pyopenms
   :sync: pyopenms

   ```{code-block} bash
   conda install pyopenms
   ```
   :::

   :::{tab-item} openms-thirdparty
   :sync: openms-thirdparty

   ```{code-block} bash
   conda install openms-thirdparty
   ```
   :::

   ::::

   or our own `openms` channel for nightly snapshots (which are build based on the same bioconda dependencies)

   ::::{tab-set}

   :::{tab-item} openms
   :sync: openms

   ```{code-block} bash 
   conda install -c openms openms
   ```
   :::

   :::{tab-item} libopenms
   :sync: libopenms

   ```{code-block} bash
   conda install -c openms libopenms
   ```
   :::

   :::{tab-item} pyopenms
   :sync: pyopenms

   ```{code-block} bash
   conda install -c openms pyopenms
   ```
   :::

   :::{tab-item} openms-thirdparty
   :sync: openms-thirdparty

   ```{code-block} bash
   conda install -c openms  openms-thirdparty
   ```
   :::

   ::::