KNIME
=====

Using OpenMS in combination with KNIME, you can create, edit, open, save, and run workflows that combine TOPP tools with
the powerful data analysis capabilities of KNIME. Workflows can be created conveniently in a graphical user interface.
The parameters of all involved tools can be edited within the application and are also saved as part of the workflow.
Furthermore, KNIME interactively performs validity checks during the workflow editing process, to make it more
difficult to create an invalid workflow. Throughout most parts of this tutorial, you will use KNIME to create and execute
workflows. The first step is to become familiar with KNIME. Additional information on the basic usage of KNIME can be
found on the KNIME [Getting Started page](https://www.knime.com/getting-started-guide). However, the most important
concepts will also be reviewed in this tutorial.

## General Remarks

- This handout will guide you through an introductory tutorial for the OpenMS/TOPP software package[^1].

- OpenMS[^2]<sup>,</sup>[^3] is a versatile open-source library for mass spectrometry data analysis. Based on this library, we offer a collection of command-line tools ready to be used by end users. These so-called TOPP tools (short for ”The OpenMS Pipeline”)[^4] can be understood as small building blocks of arbitrarily complex data analysis workflows.

- In order to facilitate workflow construction, OpenMS was integrated into KNIME[^5], the Konstanz Information Miner, an open-source integration platform providing a powerful and flexible workflow system combined with advanced data analytics, visualization, and report capabilities. Raw MS data as well as the results of data processing using TOPP can be visualized using TOPPView[^6].

- This tutorial was designed for use in a hands-on tutorial session but can also be worked through at home using the online resources. You will become familiar with some of the basic functionalities of OpenMS/TOPP, TOPPView, as well as KNIME and learn how to use a selection of TOPP tools used in the tutorial workflows.

- If you are attending the tutorial and received a USB stick, all sample data referenced in this tutorial can be found in the {path}`C:,Example_Data` folder, on the USB stick, or released online on our [Archive](https://abibuilder.cs.uni-tuebingen.de/archive/openms/Tutorials/tutorial.zip).

## Getting Started

### Installation

Before we get started, we will install OpenMS with its viewer TOPPView, KNIME and the OpenMS KNIME plugin. If you take part in a live training session you will have likely received an USB stick from us that contains the required data and software. If we provide laptops with the software you may of course skip the installation process and continue reading the next section.
If you are doing this tutorial online, choose online in the following tab(s).

:::::{tab-set}

::::{tab-item} Online
:sync: online

If you are working through this tutorial at home/online, proceed with the following steps:

- Download and install OpenMS using the installation instructions for the [OpenMS tools](/about/installation.rst).
  :::{note}
  To install the graphical application, please use the downloadable installer for your platform,
  not conda, nor docker.
  :::

- Download and install [KNIME](https://www.knime.org/downloads/)

::::

::::{tab-item} USB Stick
:sync: usb

Please choose the directory that matches your operating system and execute the installer.

For Windows, you run:

```{note}
The OpenMS installer for windows now supports installing only for a single user. If you choose this option the location of the tools will be different than {{ '{path}'+'`C:,Program Files,OpenMS-{0}`'.format(version) }} specified in this document. In most cases they will install to {{ '{path}'+'`C:,Users,$YOUR_USER,AppData,Local,OpenMS-{0}`'.format(version) }} where $YOUR_USER is replaced with your username.
```
- The OpenMS installer: {{ '{path}'+'`Windows,OpenMS-{0}-Win64.exe`'.format(version) }}
- The KNIME installer: {{ '{path}'+'`Windows,KNIME-{0}-Installer-64bit.exe`'.format(knime_version) }}

On macOS(x86), you run:

- The OpenMS installer: {{ '{path}'+'`Mac,OpenMS-{0}-macOS.dmg`'.format(version) }}
- The KNIME installer: {{ '{path}'+'`Mac,knime_{0}.app.macosx.cocoa.x86_64.dmg`'.format(knime_version) }}

On macOS(arm), you run:

- The OpenMS installer: {{ '{path}'+'`Mac,OpenMS-{0}-macOS.dmg`'.format(version) }}
- The KNIME installer: {{ '{path}'+'`Mac,knime_{0}.app.macosx.cocoa.aarch64.dmg`'.format(knime_version) }}

On Linux:
- The OpenMS package: {{ '{path}'+'`Linux,OpenMS-{0}-Debian-Linux-x86_64.deb`'.format(version) }} can be installed with your package manager
- The KNIME package can be extracted to a folder of your choice from {{ '{path}'+'`knime_{0}.linux.gtk.x86_64.tar`'.format(knime_version) }}

 ```{note}
 You can also install OpenMS via your package manager (version availability not guaranteed) or build it on your own with our
 [build instructions](/about/installation/installation-on-gnu-linux.md#build-openms-from-source).
 ```

::::
:::::

### KNIME Modern and Classic UI ####

Since version 5.0 KNIME has a new updated user interface. For the purposes of this tutorial we will continue to use the "classic user interface".
Depending on your OS KNIME may have started automatically in the Modern UI, which looks like the following:
|![ms2 spectrum](/_images/openms-user-tutorial/introduction/KNIME_switch_to_classic.png)|
|:--:|
|Figure 5.5: The modern KNIME UI. To switch back to the classic UI, select "Menu" and click "Switch to classic user interface"|

### Plugin and dependency

Before we can start with the tutorial, we need to install all the required extensions for KNIME.
Since KNIME 3.2.1, the program automatically
detects missing plugins when you open a workflow but to make sure that the right source for the
OpenMS plugin is chosen, please follow the instructions here.

#### Required KNIME plugins

First, we install some additional extensions that are
required by our OpenMS nodes or used in the Tutorials for downstream processing, visualization or reporting.

1. In KNIME, click on **Help** > **Install New Software**.
2. From the '**Work with**:' drop-down list, select the _update site_ 'KNIME 5.2 - https://update.knime.com/analytics-platform/5.2'
3. Now select the following KNIME core plugins from the KNIME & Extensions category
- KNIME Base Chemistry Types & Nodes
- KNIME Chemistry Add-Ons
- KNIME Interactive R Statistics Integration
- KNIME Report Designer
4. Click on **Next** and follow the instructions (it's not necessary to restart KNIME now).
5. Click again on **Help** > **Install New Software**
6. From the '**Work with**:' drop-down list, select the _update site_ 'KNIME Community Extensions (Trusted) - https://update.knime.com/community-contributions/trusted/5.2'
7. From the "KNIME Community Contributions - Cheminformatics" category select
- RDKit Nodes Feature
8. From the "KNIME Community Extensions - Other" category select
- Generic Worfkflow Nodes for KNIME
9.  Click on **Next**  and follow the instructions and after a restart of KNIME the dependencies will be installed.

#### R programming language and its KNIME integration

In addition, we need to install `R` for the statistical downstream analysis. Choose the directory that matches your
operating system, double-click the `R` installer and follow the instructions. We recommend to use the default settings
whenever possible. On macOS you also need to install `XQuartz` from the same directory.

Afterwards open your `R` installation. If you use Windows, you will find an ”R x64 4.3.2” icon on your desktop. If you use
macOS, you will find R in your Applications folder. In `R`, type the following lines (you might also copy them from the file
{path}`R,install_R_packages.R` on the USB stick):

```r
install.packages('Rserve',,"http://rforge.net/",type="source")
install.packages("Cairo")

install.packages("devtools")
install.packages("ggplot2")
install.packages("ggfortify")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install()
BiocManager::install(c("MSstats"))
```
In KNIME, click on **File** > **Preferences**, select the category **KNIME** > **R** and set the ”Path to R Home” to
your installation path. You can use the following settings, if you installed R as described above:

- Windows: `C:\Program Files\R\R-4.3.2'`
- macOS: `/Library/Frameworks/R.framework/Versions/4.3/Resources`

#### KNIME OpenMS plugin

You are now ready to install the OpenMS nodes.
- In KNIME, click on **Help** > **Install New Software**

You now have to choose an _update site_ to install the OpenMS plugin from. Which _update site_ to choose depends on if you received an USB stick
in a hands-on Tutorial or if you are doing this Tutorial online.

::::{tab-set}

:::{tab-item} Online
:sync: online

To install the OpenMS KNIME plugin from the internet, do the following:

1. From the '**Work with**:' drop-down list, select the _update site_ 'KNIME Community Extensions (Trusted) - https://update.knime.com/community-contributions/trusted/5.2'
2. Now select the following plugin from the "KNIME Community Contributions - Bioinformatics & NGS" category
- OpenMS
- OpenMSThirdParty
3. Click on **Next**  and follow the instructions and after a restart of KNIME the OpenMS nodes will be available in the Node repository under
   "Community Nodes".

```{note}
If this does not work for you, report it and while waiting for a reply/fix, try to use an _update site_ of an older KNIME version by editing the KNIME
version number in the URL or by using our inofficial _update site_ at https://abibuilder.cs.uni-tuebingen.de/archive/openms/knime-plugin/updateSite/release/latest
```

:::

:::{tab-item} USB
:sync: usb

We included a custom KNIME update site to install the OpenMS KNIME plugins from the USB stick. If you do not have a stick available, please see below.

- In the now open dialog choose **Add** (in the upper right corner of the dialog) to define a new update site. In the
  opening dialog enter the following details.

  **Name:** OpenMS {{ version }} UpdateSite

  **Location:** {{ '`file:/KNIMEUpdateSite/{0}/`'.format(version) }}
- After pressing **OK** KNIME will show you all the contents of the added Update Site.

```{note}
From now on, you can use this repository for plugins in the **Work with**: drop-down list.
```

- Select the OpenMS nodes in the ”Uncategorized” category and click **Next**.
- Follow the instructions and after a restart of KNIME the OpenMS nodes will be available in the Node repository under
  "Community Nodes”.

:::

:::{tab-item} Online experimental
:sync: onlineexp

To install the nightly/experimental version of the OpenMS KNIME plugin from the internet, do the following:

- In the now open dialog, choose **Add** (in the upper right corner of the dialog) to define a new _update site_. In the opening dialog enter the following details.

  **Name:** OpenMS {{ version }} UpdateSite

  **Location:** https://abibuilder.cs.uni-tuebingen.de/archive/openms/knime-plugin/updateSite/nightly/
- After pressing **OK** KNIME will show you all the contents of the added Update Site.

```{note}
From now on, you can use this repository for plugins in the **Work with:** drop-drown list.
```
- Select the OpenMS nodes in the "Uncategorized" category and click **Next**.
- Follow the instructions and after a restart of KNIME the OpenMS nodes will be available in the Node repository under
  "Community Nodes".

:::
::::


```{toctree}
---
maxdepth: 1
titlesonly: true
---

knime-user-tutorial/file-conversion.md
knime-user-tutorial/knime-gui.md
knime-user-tutorial/minimal-workflow.md
knime-user-tutorial/lfq-peptide-protein.md
knime-user-tutorial/msstats.md
knime-user-tutorial/lfq-metabolites.md
knime-user-tutorial/openswath.md
knime-user-tutorial/openswath-metabolomics.md
knime-user-tutorial/quality-control.md

```

## References

[^1]: OpenMS, <a href="http://www.openms.de/">OpenMS home page</a> [online].

[^2]: M. Sturm, A. Bertsch, C. Gröpl, A. Hildebrandt, R. Hussong, E. Lange, N. Pfeifer,
O. Schulz-Trieglaff, A. Zerck, K. Reinert, and O. Kohlbacher, <a href="http://dx.doi.org/10.1186/1471-2105-9-163">OpenMS - an opensource software framework for mass spectrometry</a>., BMC bioinformatics 9(1)
(2008), <a href="https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-163">doi:10.1186/1471-2105-9-163</a>. 7, 83

[^3]: H. L. Röst, T. Sachsenberg, S. Aiche, C. Bielow, H. Weisser, F. Aicheler, S. Andreotti,
H.-C. Ehrlich, P. Gutenbrunner, E. Kenar, et al., OpenMS: a flexible open-source
software platform for mass spectrometry data analysis, Nature Methods 13(9),
741–748 (2016). 7

[^4]: O. Kohlbacher, K. Reinert, C. Gröpl, E. Lange, N. Pfeifer, O. Schulz-Trieglaff, and
M. Sturm, <a href="http://view.ncbi.nlm.nih.gov/pubmed/17237091">TOPP–the OpenMS proteomics pipeline</a>., Bioinformatics 23(2) (Jan.
2007). 7, 83

[^5]: M. R. Berthold, N. Cebron, F. Dill, T. R. Gabriel, T. Kötter, T. Meinl, P. Ohl, C. Sieb, K. Thiel, and B. Wiswedel, KNIME: The Konstanz Information Miner, in Studies in Classification, Data Analysis, and Knowledge Organization (GfKL 2007), Springer, 2007.

[^6]: M. Sturm and O. Kohlbacher, <a href="http://dx.doi.org/10.1021/pr900171m">TOPPView: An Open-Source Viewer for Mass Spectrometry Data</a>, Journal of proteome research 8(7), 3760–3763 (July 2009), <a href="https://pubs.acs.org/doi/10.1021/pr900171m?cookieSet=1">doi:10.1021/pr900171m</a>. 7
