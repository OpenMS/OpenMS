TOPPAS
======

TOPPAS is the build-in workflow editor of OpenMS.
All TOPP tools can be chained, configured and executed.

TOPPAS workflows run on the local machine where TOPPAS is executed and thus only scale according to the hardware at hand. No automatic distribution across a cluster is supported.
TOPPAS runs OpenMS TOPP tools only, including the adapters that wrap search engines and other third-party programs;
other external tools cannot be called from a TOPPAS workflow. Hand the resulting data to such tools, or use a workflow
system such as Nextflow or Galaxy, which can integrate any tool.

The strong point of TOPPAS is that it ships with OpenMS natively in the graphical and command-line tools.
It also has a very shallow learning curve, making it very intuitive to create workflows.

See the [TOPPAS tutorial](https://openms.de/doxygen/nightly/html/TOPPAS_tutorial.html) for more details.