Contribute
==========

## Reporting Bugs and Issues

A list of known issues in the current OpenMS release can be found [here](https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/nightly/html/known_dev_bugs.html). 
Please check if your OpenMS version matches the current version and if the bug has already been reported.

In order to report a new bug, please create a [GitHub issue](manual/contribute.md#Write and Label GitHub Issues) or [contact us](/about/communication.md).

Include the following information in your bug report:

1. The command line (i.e. call) including the TOPP tool and the arguments you used, or the steps you followed in a GUI
   tool (e.g. TOPPView) - e.g. `FeatureFinderCentroided -in myfile.mzML -out myfile.featureXML`.
2. The output of OpenMS/TOPP (or a screenshot in case of a GUI problem).
3. Operating system (e.g. "Windows XP 32 bit", "Win 7 64 bit", "Fedora 8 32 bit", "macOS 10.6 64 bit").
4. OpenMS version (e.g. "OpenMS 1.11.1", "Revision 63082 from the SVN repository").
5. OpenMS architecture ("32 bit" or "64 bit")

Please provide files that we need to reproduce the bug (e.g. `TOPP INI` files, data files — usually mzML) via a download
link, via the mailing list or by directly contacting one of the developers.


## Write and Label GitHub Issues

### Create an Issue

To create an issue:

1. Go to the [OpenMS codebase](https://github.com/OpenMS/OpenMS).
2. Submit an [issue](https://github.com/OpenMS/OpenMS/issues/new).

The issue will be listed under **Issues**.

### Label an Issue

To label an issue:

1. On the right of the screen, select the cog icon under **Labels**.
2. Choose a label from the list. Normally, an issue can have one or more of the following labels:
   - **defect**: A defect refers to a bug in OpenMS. This is a high priority issue.
   - **enhancement**: An enhancement refers to a feature idea to enhance the current OpenMS code. This is a medium
     priority issue.
   - **task**: A task refers to a single piece of work that a developer can undertake. This is a medium priority issue.
   - **refactoring**: A refactoring issue refers to a suggestion to streamline the code without changing how the code
     function.
   - **question**: A question could trigger to a discussion about tools, parameters and scientific tasks.



```{toctree}
:maxdepth: 1

contribute/pull-request-checklist.md
contribute/openms-git-workflow.md

```