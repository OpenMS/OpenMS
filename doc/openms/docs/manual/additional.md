Additional
==========

## Graphical User Interfaces

OpenMS provides additional graphical user interfaces besides TOPPAS and TOPPView, designed for users who want easy access to the OpenMS tools. These interfaces include:

- **INIFileEditor**

  A GUI application used to edit INI files. INI files are used to configure the parameters of OpenMS tools. INI files are files with the extension `.ini`. For mor information, read our [INIFileEditor](additional/ini-file-editor.md) section.


A possible workflow would consist of the following steps:

1. Generate an INI file from the [command line](/getting-started/topp-tools.md#command-line-interface).
2. Edit the INI file in the INIFile Editor.
3. Import data into TOPPView.
4. Apply TOPP tool to data in TOPPView. You will need to load the INI file edited in step 1.


```{toctree}
---
maxdepth: 1
hidden: True
---

additional/ini-file-editor.md
```