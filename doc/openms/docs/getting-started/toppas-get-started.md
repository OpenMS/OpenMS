TOPPAS
======

TOPPAS is the build-in workflow editor of OpenMS.
All TOPP tools can be chained, configured and executed.

TOPPAS workflows run on the local machine where TOPPAS is executed and thus only scale according to the hardware at hand. No automatic distribution across a cluster is supported.
TOPPAS runs OpenMS TOPP tools only, including the adapters that wrap search engines and other third-party programs;
other external tools cannot be called from a TOPPAS workflow. Hand the resulting data to such tools, or use a workflow
system such as Nextflow or Galaxy, which can integrate any tool.

The strong point of TOPPAS is that it is made for OpenMS: it runs the TOPP tools of the OpenMS installation it is
installed into, and the same workflows run on the command line with ExecutePipeline.
It also has a very shallow learning curve, making it very intuitive to create workflows.

TOPPAS is part of the graphical tools of OpenMS, which are developed in
[OpenMS-TOPP-GUI](https://github.com/OpenMS/OpenMS-TOPP-GUI) and installed on top of OpenMS
(see [Installation](/about/installation.rst)).

See the [TOPPAS tutorial](https://archive.openms.de/openms/Documentation/release/latest/html-gui/TOPPAS_tutorial.html) for more details.