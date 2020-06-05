#************************************ Enter script definitions **********************************#

!define APP_NAME "Advanced Uninstall Log Demo"
#
#
;..................................................................................................
;Following two definitions required. Uninstall log will use these definitions.
;You may use these definitions also, when you want to set up the InstallDirRagKey,
;store the language selection, store Start Menu folder etc.
;Enter the windows uninstall reg sub key to add uninstall information to Add/Remove Programs also.

!define INSTDIR_REG_ROOT "HKLM"
!define INSTDIR_REG_KEY "Software\Microsoft\Windows\CurrentVersion\Uninstall\${APP_NAME}"
;..................................................................................................

#************************************* Include NSIS Headers ***********************************#
#
#
;..................................................................................................
;include the Uninstall log header
!include AdvUninstLog.nsh
;..................................................................................................

Name "${APP_NAME}"
OutFile "${APP_NAME}[DEFAULT UI].exe"
ShowInstDetails show
ShowUninstDetails show
LicenseData "${NSISDIR}\License.txt"
LicenseBkColor 0xFFFFFF
InstallDir "$PROGRAMFILES\${APP_NAME}"
InstallDirRegKey ${INSTDIR_REG_ROOT} "${INSTDIR_REG_KEY}" "InstallDir"

;..................................................................................................
;Specify the preferred uninstaller operation mode, either unattended or interactive.
;You have to type either !insertmacro UNATTENDED_UNINSTALL, or !insertmacro INTERACTIVE_UNINSTALL.
;Be aware only one of the following two macros has to be inserted, neither both, neither none.

  !insertmacro UNATTENDED_UNINSTALL
  ;!insertmacro INTERACTIVE_UNINSTALL
;..................................................................................................

Page License
Page Components
Page Directory 
Page InstFiles


Section "Main Application" sec01

        SetOutPath '$INSTDIR'

;After set the output path open the uninstall log macros block and add files/dirs with File /r
;This should be repeated every time the parent output path is changed either within the same
;section, or if there are more sections including optional components.
        !insertmacro UNINSTALL.LOG_OPEN_INSTALL

        File /r /x "uninst-nsis.exe" /x "Docs" /x "Stubs" "${NSISDIR}\*"

;........................................................................................
;Now for the 2nd run (update installation mode) uncomment the File /r instruction below,
;and comment the File /r instruction above. New installation files would be appended to log.
;Everything that's included within the macros block uninstall_log_open/close added in log.
;Files that added outside the block either manually or from application's usage,
;never added to log, so uninstaller knows about them and requests confirmation to delete.

        ;File /r '${NSISDIR}\Docs'
;........................................................................................

;Once required files/dirs added and before change the parent output directory we need to
;close the opened previously uninstall log macros block.
        !insertmacro UNINSTALL.LOG_CLOSE_INSTALL

        CreateDirectory '$SMPROGRAMS\${APP_NAME}'
        CreateShortcut '$SMPROGRAMS\${APP_NAME}\nsis.lnk' '$INSTDIR\nsis.exe'
        ;create shortcut for uninstaller always use ${UNINST_EXE} instead of uninstall.exe
        CreateShortcut '$SMPROGRAMS\${APP_NAME}\uninstall.lnk' '${UNINST_EXE}'

        WriteRegStr ${INSTDIR_REG_ROOT} "${INSTDIR_REG_KEY}" "InstallDir" "$INSTDIR"
        WriteRegStr ${INSTDIR_REG_ROOT} "${INSTDIR_REG_KEY}" "DisplayName" "${APP_NAME}"
        ;Same as create shortcut you need to use ${UNINST_EXE} instead of anything else.
        WriteRegStr ${INSTDIR_REG_ROOT} "${INSTDIR_REG_KEY}" "UninstallString" "${UNINST_EXE}"

#...............................................................................
;instead of adding files/folders manually to perform our test, we do it here :)
        File /r '${NSISDIR}\Stubs'
#...............................................................................

SectionEnd


Section "Application Data" sec02

        SetOutPath "$APPDATA\${APP_NAME}"

;After set the output path open the uninstall log macros block and add files/dirs with File /r
;This should be repeated every time the parent output path is changed either within the same 
;section, or if there are more sections including optional components.
        !insertmacro UNINSTALL.LOG_OPEN_INSTALL

        File /r '${NSISDIR}\Stubs\*'

;........................................................................................
;Now for the 2nd run (update installation mode) uncomment the File /r instruction below,
;and comment the File /r instruction above. New installation files would be appended to log.
;Everything that's included within the macros block uninstall_log_open/close added in log.
;Files that added outside the block either manually or from application's usage,
;never added to log, so uninstaller knows about them and requests confirmation to delete.

        ;File /r '${NSISDIR}\Contrib'
;........................................................................................

;Once required files/dirs added and before change the parent output directory we need to
;close the opened previously uninstall log macros block.
        !insertmacro UNINSTALL.LOG_CLOSE_INSTALL

#...............................................................................
;instead of adding files/folders manually to perform our test, we do it here :)
        File /r '${NSISDIR}\Plugins'
#...............................................................................

SectionEnd


Function .onInit

        ;prepare log always within .onInit function
        !insertmacro UNINSTALL.LOG_PREPARE_INSTALL

FunctionEnd


Function .onInstSuccess

         ;create/update log always within .onInstSuccess function
         !insertmacro UNINSTALL.LOG_UPDATE_INSTALL

FunctionEnd

#######################################################################################

Section UnInstall

         ;begin uninstall, especially for MUI could be added in UN.onInit function instead
         !insertmacro UNINSTALL.LOG_BEGIN_UNINSTALL

         ;uninstall from path, must be repeated for every install logged path individual
         !insertmacro UNINSTALL.LOG_UNINSTALL "$INSTDIR"

         ;uninstall from path, must be repeated for every install logged path individual
         !insertmacro UNINSTALL.LOG_UNINSTALL "$APPDATA\${APP_NAME}"

         ;end uninstall, after uninstall from all logged paths has been performed
         !insertmacro UNINSTALL.LOG_END_UNINSTALL

        Delete "$SMPROGRAMS\${APP_NAME}\nsis.lnk"
        Delete "$SMPROGRAMS\${APP_NAME}\uninstall.lnk"
        RmDir "$SMPROGRAMS\${APP_NAME}"

        DeleteRegKey /ifempty ${INSTDIR_REG_ROOT} "${INSTDIR_REG_KEY}"

SectionEnd

/*
Function UN.onInit

         ;begin uninstall, could be added on top of uninstall section instead
         !insertmacro UNINSTALL.LOG_BEGIN_UNINSTALL

FunctionEnd
*/
