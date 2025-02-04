
##############################
### register filetypes (code borrowed from VLC project)
;; Function that register one extension for VLC

Function RegisterExtension
  ; back up old value for extension $R0 (eg. ".mzData") and Program $R1 (e.g. TOPPView)
  ReadRegStr $1 HKCR "$R0" ""
  StrCmp $1 "" NoBackup
    StrCmp $1 "OPENMS$R0" "NoBackup"
    WriteRegStr HKCR "$R0" "OPENMS.backup" $1
NoBackup:
  WriteRegStr HKCR "$R0" "" "OPENMS$R0"
  ReadRegStr $0 HKCR "OPENMS$R0" ""
  WriteRegStr HKCR "OPENMS$R0" "" "OPENMS data file ($R0)"
  #WriteRegStr HKCR "OPENMS$R0\shell" "" "Open"
  WriteRegStr HKCR "OPENMS$R0\shell\open\command" "" '"$INSTDIR\bin\$R1.exe" "%1"'
  WriteRegStr HKCR "OPENMS$R0\DefaultIcon" "" '"$INSTDIR\share\OpenMS_$R1.ico",0'

;;; Vista Only part
  ; Vista detection
#  ReadRegStr $R1 HKLM "SOFTWARE\Microsoft\Windows NT\CurrentVersion" CurrentVersion
#  StrCpy $R2 $R1 3
#  StrCmp $R2 '6.0' ForVista ToEnd
#ForVista:
#  WriteRegStr HKLM "Software\Clients\Media\VLC\Capabilities\FileAssociations" "$R0" "VLC$R0"

ToEnd:
FunctionEnd

;; Function that removes one extension that OPENMS owns.
Function un.RegisterExtension
  ;start of restore script
  ReadRegStr $1 HKCR "$R0" ""
  StrCmp $1 "OPENMS$R0" 0 NoOwn ; only do this if we own it
    ; Read the old value from Backup
    ReadRegStr $1 HKCR "$R0" "OPENMS.backup"
    StrCmp $1 "" 0 Restore ; if backup="" then delete the whole key
      DeleteRegKey HKCR "$R0"
    Goto NoOwn
Restore:
      WriteRegStr HKCR "$R0" "" $1
      DeleteRegValue HKCR "$R0" "OPENMS.backup"
NoOwn:
    DeleteRegKey HKCR "OPENMS$R0" ;Delete key with association settings
#   DeleteRegKey HKLM "Software\Clients\Media\OPENMS\Capabilities\FileAssociations\OPENMS$R0" ; for vista
FunctionEnd

!macro RegisterExtensionSection EXT PROGRAMEXE
  Section ${EXT}
    SectionIn 1 3
    Push $R0
		Push $R1
		StrCpy $R0 ${EXT}
		StrCpy $R1 ${PROGRAMEXE}
		Call RegisterExtension
		Pop $R1
		Pop $R0
  SectionEnd
!macroend

!macro UnRegisterExtensionSection EXT PROGRAMEXE
  Push $R0
	Push $R1
  StrCpy $R0 ${EXT}
	StrCpy $R1 ${PROGRAMEXE}
  Call un.RegisterExtension
	Pop $R1
  Pop $R0
!macroend


!macro REBOOT_ON_INCOMPLETE_DELETION PROGRAMPATH
    ClearErrors
    SetRebootFlag false
    # old installer directory structure
    Delete /REBOOTOK ${PROGRAMPATH}\lib\*.dll
    Delete /REBOOTOK ${PROGRAMPATH}\GUI\*.exe
    Delete /REBOOTOK ${PROGRAMPATH}\TOPP\*.exe
    # new installer directory structure
    Delete /REBOOTOK ${PROGRAMPATH}\bin\*.dll
    Delete /REBOOTOK ${PROGRAMPATH}\bin\*.exe

    IfRebootFlag 0 noreboot
        MessageBox MB_OK|MB_ICONQUESTION "Not all OpenMS files could be deleted. A reboot is required! \
                                             Click 'OK' when you are ready to reboot."
        Reboot
    noreboot:
!macroend

