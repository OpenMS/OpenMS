Adding New Tool to The TOPP suite
=====================================

## The OpenMS pipeline (TOPP)

Any tool that is written with the OpenMS library can easily be made into a TOPP tool by simply using the OpenMS command
line parser which is able to parse ParamXML, a powerful XML based description of the tool. Hence most analysis algorithms
in OpenMS are available as a stand-alone tool which can be called on the command line or integrated into workflow engines
via the CTD mechanism. A current list of TOPP tools can be found in [the documentation](https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/TOPP_documentation.html).

## What do I have to do to add a new TOPP tool?

The recommended way is to inherit from the class TOPPBase as in existing TOPP tools (sources available in /src/topp/). This will add command line parsing functionality to your tool as described in the TOPP section of this page.

- Add the code to `src/topp/` and register it in `src/topp/executables.cmake`
- Add your tool (with the correct category) to `getTOPPToolList()` in `src/openms/source/APPLICATIONS/ToolHandler.cpp`.
  This creates a doxygen page with the `–help` output of the tool (using `TOPPDocumenter`). This page must be included
  at the end of the doxygen documentation of your tool (see other tools for an example).
- Add it to the TOPP docu page (in `doc/doxygen/public/TOPP.doxygen`)
- Add the name to `src/topp/executables.cmake`
- Write a TOPP test (add it to `src/tests/topp/CMakeLists.txt`)

```{warning}
Handle any kind of input files to your TOPP tool via command line flags and use the `${DATA_DIR_TOPP}` prefix. Use
ini-files to specify output-files, but not input-files. Doing otherwise will break out-of-source builds.
```

```{hint}
add `-test` to the call of your TOPP tool and also create the expected output that you put in `src/tests/topp` with that
flag active. The flag ensures that UniqueId's,  dates etc are equal no matter where and when the tool is run.
```

## I want to implement a new file adapter. What is to be done?

First, add a file adapter class to the `include/OpenMS/FORMAT/` and `source/FORMAT/` folders. The file adapter should
implement a default constructor, a load method and a store method. Make sure your code conforms to the OpenMS Coding
conventions. For automatic file type recognition, you need to

- register your new file type at the Type enum in `/include/OpenMS/FORMAT/FileTypes.h`,
- flag the file type as supported in the isSupported method of `/source/FORMAT/FileHandler.C`
- register the file extension in the getTypeByFileName method of `/source/FORMAT/FileHandler.C`

If the new file is a peak or feature file format you should also add it to loadExperiment or loadFeatures, respectively,
of the FileHandler class. To add the file format to the TOPPView open dialog, you have to modify the file
`/source/APPLICATIONS/TOPPViewBase.C`.

- Add the file extensions to the filter_all and filter_single variables of the getFileList_ method.

To add your format to TOPP applications:

- add the file extension to the extensions list of the respective parameter:
  ```
  e.g. setValidStrings_("in_type", StringList::create("mzData,mzXML,mzML")); in FileInfo
  ```

## How to create an icon file for a TOPP tool under Windows?

- Create an .ico file: first, you need some graphics program (The GIMP is recommended) think of a motive and remind
  yourself that you have limited space. Create at least a 16x16, 32x32, 48x48 and 64x64 pixel version and save each of
  them in a separate layer of the respective size. Do not add any larger sized layers, since Win XP will not display any
  icon then. When saving the image as type `.ico` the GIMP will ask you for the color depth of each layer. As it is
  recommended to have multiple color depths of each icon-size, go back to the layers and duplicate each layer twice.
  That should give you 12 layers. Now, save the image as `.ico` (e.g. TOPPView.ico) file, giving each group of equal
  sized layers a 32 bit (8 bit transparency), 8 bit (1 bit transparency), 4 bit (1 bit transparency) color depth.

```{attention}
Make sure to assign the higher color depth to the upper layers as Windows will not pick the highest possible color
otherwise.
```

- Create a resource file: Create a text file named `.rc` (e.g. TOPPView.rc) Insert the following line: 101 ICON
  "TOPPView.ico" , replacing TOPPView with your binary name. Put both files in `OpenMS/source/APPLICATIONS/TOPP/`
  (similar files for other TOPP tools already present). Re-run cmake and re-link your TOPP tool. 

Voila. You should have an iconized TOPP tool.

## Develop your Tool in an external project using OpenMS

To include the OpenMS library in one of your projects, we recommend to have a look at a small emulated external project
in our repository. We strongly suggest to use CMake for building your project together with OpenMS to  make use of the
macros and environment information generated during the build of the OpenMS library.

## The Common Tool Description (CTD)

The CTD is a format developed from the OpenMS team to allow the user to use TOPP tools also in other workflow engines.
Each tool can output a CTD description of itself (the XML scheme for the CTD can be found here), which can then be used
by a node generator program to generate nodes for different workflow engines. The CTD mechanism is shared by OpenMS with
other mature libraries like SeqAn and BALL. An example for a node generation program are the Generic KNIME Nodes. The
most complete description on how to generate your own Generic KNIME Nodes based on a CTD (e.g. from your freshly
developed command line tool), can be found on the SeqAn documentation. We are working on a tutorial specifically
tailored to OpenMS.
