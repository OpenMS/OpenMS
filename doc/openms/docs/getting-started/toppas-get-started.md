TOPPAS
======

TOPPAS is the build-in workflow editor of OpenMS.
All TOPP tools can be chained, configured and executed.

TOPPAS workflows run on the local machine where TOPPAS is executed and thus only scale according to the hardware at hand. No automatic distribution across a cluster is supported.
Also, external tools (anything other than OpenMS TOPP tools), can only be called using a special 'GenericWrapper' node. We generally recommend to run only OpenMS-specific tools in TOPPAS
and hand the resulting data to other tools, or use other workflow systems, such as KNIME which can fully integrate other tools.

The strong point of TOPPAS is that it ships with OpenMS natively in the graphical and command-line tools.
It also has a very shallow learning curve, making it very intuitive to create workflows.

See the [TOPPAS tutorial](https://openms.de/doxygen/nightly/html/TOPPAS_tutorial.html) for more details.