This ReadMe describes how OpenMS deals with TOPP tool description (*.ttd) files and the external tools wrapped by them.

== General ==

Using *.ttd files OpenMS can be instructed to run external tools (like msconvert from Proteowizard) in a convenient way.
Each *.ttd file can contain one or more <tool> sections, each telling OpenMS how to call the respective tool.
You can then use our TOPP tool

  GenericWrapper

to run the external tool.
GenericWrapper will dynamically change its available '-type' parameter list, depending on the *.ttd files present in this folder.

The name of the *.ttd file does not really matter, but it should be descriptive.

This way you can now easily integrate external tools into TOPPAS, by just adding a GenericWrapper node of the desired '-type'.

===============
== For Users ==
===============

If you want to add a custom *.ttd file (obtainable from our website or other users), simply add it to this folder.
Some tools (like 'mail' on linux), might only be available on certain Operating Systems. In this case you should
place them in the respective subfolder. The folder of this ReadMe and the subfolder matching your Operating System will be
scanned for *.ttd files.
Additionally you can create an environment variable OPENMS_TTD_PATH and use it to point to one (!) custom directory of your choice.
This is useful when upgrading OpenMS to another version while keeping all your custom *.ttd files.

If you encounter a bug, file a bug report via https://github.com/OpenMS/OpenMS/issues
and provide the output of GenericWrapper (ideally with '-debug 10' option enabled)
and your input files (if applicable).

After adding/removing *.ttd files you need to restart any open instances of TOPPAS to see the updated list of types for GenericWrapper.


===========================
== For experienced Users ==
===========================

Hints -- read(!):
 - When mapping a SINGLE filename while building the command line in 'cloptions', wrap them in QUOTES! (to avoid the 'spaces in filenames' problem)
 - When mapping a LIST of filenames (from ITEMLIST parameters) while building the command line in 'cloptions', NO NOT QUOTE! (GenericWrapper will do this automatically)
 - Find a good name for <type> - it should be descriptive duplicate names are not allowed

If you want to wrap your own tools and write custom *.ttd files,
a good starting point is the 

  TEMPLATE.ttd_

file. See description inside the file.

The template's suffix is 'ttd_' in order to disregard it as an active *.ttd file.
Make a copy of it and modify to your needs.

Once you are done, place the new *.ttd in either this folder (for all Operating Systems) or the correct subfolder.

You can also have a look at existing wrappers and get inspiration.

A note on using more complicated shell commands (e.g. piping):
  - as we use QProcess to call the external tool, no piping and shell internal commands (such as 'echo' on Windows) are supported
    To work around this you have to specify the command tool itself, e.g.
	
	<path>cmd</path>
	<cloptions>/C "echo ""bla bla"" > testfile"</cloptions>
	
	or
	
	<path>sh</path>
	<cloptions>-c "echo """%2""" > testfile"</cloptions>
                  
	Note that escaping of internal quotes works different on each platform. You need to test this!
	
