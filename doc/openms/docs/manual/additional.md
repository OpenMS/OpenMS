Additional
==========

## Graphical User Interfaces

OpenMS provides additional graphical user interfaces besides TOPPAS and TOPPView, designed for users who want easy access to TOPP tools. These interfaces include:

- **INIFileEditor**

  A GUI application used to edit TOPP INI files. TOPP INI files are used to configure TOPP tool parameters. TOPP INI files are files with the extension `.ini`. For mor information, read our [INIFileEditor](additional/ini-file-editor.md) section.

- **SwathWizard**
  An application for SWATH analysis. SwathWizard is used to analyze DIA swath data. For more information, read our [SwathWizard](additional/swathwizard.md) section.


A possible workflow would consist of the following steps:

1. Generate a TOPP INI file from the [command line](/getting-started/topp-tools.md#command-line-interface).
2. Edit the TOPP INI file in the INIFile Editor.
3. Import data into TOPPView.
4. Apply TOPP tool to data in TOPPView. You will need to load the TOPP INI file edited in step 1.


```{toctree}
---
maxdepth: 1
hidden: True
---

additional/swathwizard.md
additional/ini-file-editor.md
```