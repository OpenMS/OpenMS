In the os_compiler_matrix.tsv please specify
the operating system and compiler combinations
that should be used to build packages/plugins
for this commit.
Only change if you know what you are doing.

Use Jenkins labels from the full matrix
(e.g. see here: https://abibuilder.informatik.uni-tuebingen.de/jenkins/job/openms_release_packaging/)
First column must be the OS, all following
columns are the compilers used on that OS.
Combinations that do not have a slave that
can run them are ignored.
Usually you want exactly one compiler for
Debian, Windows and macOS. Multiple compilers
are not yet tested.
