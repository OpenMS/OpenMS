Windows
=======

## Install via Windows installer

To Install the binary package of OpenMS & {term}`TOPP`:

1. Download the installer `OpenMS-<version>-Win64.exe` from the [archive](https://archive.openms.de/openms/OpenMSInstaller/release/latest/)
2. Execute the installer under the user account that later runs OpenMS and follow its instructions.
   
   You may see a Windows Defender Warning, since our installer is not digitally signed.
   
   Click on "More Info", and then "Run anyways".
   
   ![](/_images/installations/win/smartscreen.gif)

   When asked for an admin authentication, please enter the credentials (it is not advised to directly invoke the installer using an admin account).

```{tip}
The windows installer works with Windows 10 and 11 (older versions might still work but are untested).
```

## Known issues

1. During installation, an error message pops up, saying:
   > "The installation of the Visual Studio redistributable package ... failed. ..."

   This is a known issue with a Microsoft package, we cannot do anything about it.
   The error message will give the location where the redistributable package was extracted to. Go to this folder and
   run the executable (usually named `vcredistXXXX.exe`) as an administrator (right-click and then select **Run-As**). You will likely
   receive an error message (this is also the reason why the OpenMS setup complained about it). You might have to find
   the solution to fix the problem in your local machine. If you're lucky the error message is instructive and the
   problem is easy to fix.
2. During installation, an error message pops up saying:
   >"Error opening installation log file"

   To fix, check the system environment variables. Make sure they are apt. There should a `TMP` and a `TEMP` variable,
   and both should contain one directory only, which exists and is writable. Fix accordingly (search the internet on
   how to change environment variables on Windows).

Since OpenMS 3.6 the installer no longer includes ProteoWizard. To convert other vendor formats with `msconvert`,
install ProteoWizard from its [website](https://proteowizard.sourceforge.io/).

## Reading Thermo Fisher RAW files

On Windows, FileConverter converts `.raw` files with ThermoRawFileParser by default. The
installer puts it on the `PATH`, and it runs on the .NET Framework that Windows includes. The
other tools that read `.raw` files, and FileConverter with `-RawToMzML:reader inprocess`, use
the openms-thermo-bridge that is built into OpenMS. This requires a **.NET 8 runtime** to be
present at run time so that the managed bridge libraries can be loaded. This is the modern,
cross-platform .NET runtime and is **not** the .NET Framework that Windows includes. Which
tools read `.raw` files is described in
[Vendor formats](/getting-started/vendor-formats.md).

Download and install it from the [.NET download page](https://dotnet.microsoft.com/download).
The official installer registers the runtime globally, so no further configuration is normally
needed.

If you installed the runtime to a non-standard location (for example an xcopy install of the
.NET runtime), set the `DOTNET_ROOT` environment variable to the install directory — the folder
that contains `dotnet.exe` and the `shared\` sub-directory — so the bridge can locate the
runtime, e.g. `DOTNET_ROOT=C:\Program Files\dotnet`.
