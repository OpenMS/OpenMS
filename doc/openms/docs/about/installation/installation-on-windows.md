Windows
=======

## Install via Windows installer

To Install the binary package of OpenMS & {term}`TOPP`:

1. Download the installer `OpenMS-<version>-Win64.exe` from the [archive](https://abibuilder.cs.uni-tuebingen.de/archive/openms/OpenMSInstaller/release/latest/) 
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
   >"The installation of the Microsoft .NET 3.5 SP1' package failed!

   You must download and install it manually in order for Proteowizard to work.
   This should only happen if installion is done by selecting the "Third Party - Proteowizard" components. The reason is
   usually that **.NET 3.5 SP1** is already installed (see Windows Control Panel). If it's not installed, follow the
   instructions of the error message.
2. During installation, an error message pops up, saying:
   > "The installation of the Visual Studio redistributable package ... failed. ..."

   This is a known issue with a Microsoft package, we cannot do anything about it.
   The error message will give the location where the redistributable package was extracted to. Go to this folder and
   run the executable (usually named `vcredistXXXX.exe`) as an administrator (right-click and then select **Run-As**). You will likely
   receive an error message (this is also the reason why the OpenMS setup complained about it). You might have to find
   the solution to fix the problem in your local machine. If you're lucky the error message is instructive and the
   problem is easy to fix.
3. During installation, an error message pops up saying:
   >"Error opening installation log file"

   To fix, check the system environment variables. Make sure they are apt. There should a `TMP` and a `TEMP` variable,
   and both should contain one directory only, which exists and is writable. Fix accordingly (search the internet on
   how to change environment variables on Windows).
4. For Win8 or later, Windows will report an error while installing `.net4` as it's mostly included. But it might occur
   that `.net3.5` does not get properly installed during the process.

   Fix is to enable the .NET Framework 3.5 yourself through Control Panel. See this [Microsoft help page](https://docs.microsoft.com/en-us/dotnet/framework/install/dotnet-35-windows).aspx#ControlPanel) for detailed information. Even if this step fails, this does not affect the functionality of OpenMS, except for the executability of included third party tools (ProteoWizard).
