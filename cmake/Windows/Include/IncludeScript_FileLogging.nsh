## most of this is copied from VLC install script
## , with bugfixes and decorated with some additional comments for better understanding
!ifndef FOLDERLIST_EXE
  !define FOLDERLIST_EXE FolderList.exe
!endif

;;;;;;;;;;;;;;;
; 4. Logging  ;
;;;;;;;;;;;;;;;
Var UninstallLog   # handle to a log file
!macro OpenUninstallLog
  FileOpen $UninstallLog "$INSTDIR\uninstall.log" a
  FileSeek $UninstallLog 0 END
!macroend

!macro CloseUninstallLog
  FileClose $UninstallLog
  SetFileAttributes "$INSTDIR\uninstall.log" HIDDEN
!macroend


;;;;;;;;;;;;;;;;;;;;
; 5. Installations ;
;;;;;;;;;;;;;;;;;;;;
!macro InstallFile FILEREGEX
  # do the normal job
  File "${FILEREGEX}"
  ## add files to uninstall.log
  
  #strip source path (retain filename only)
  push $R0
  ${GetFileName} "${FILEREGEX}" $R0
  
  !define Index 'Line${__LINE__}'
  # search for filename in target path
  FindFirst $0 $1 "$OUTDIR\$R0"
  StrCmp $0 "" "${Index}-End"
  "${Index}-Loop:"
    StrCmp $1 "" "${Index}-End"
    FileWrite $UninstallLog "$OUTDIR\$1$\r$\n"
    FindNext $0 $1
    Goto "${Index}-Loop"
  "${Index}-End:"
  !undef Index
  
  pop $R0
!macroend


!macro InstallFolder FOLDER EXCLUDEREGEX
  #old: File /r "${FOLDER}"
  !system "${FOLDERLIST_EXE} /r $\"${FOLDER} $\" $\"${EXCLUDEREGEX}$\"" 
  !include "${FOLDERLIST_EXE}.txt" 
  
  Push "$OUTDIR"
  Call InstallFolderInternal
!macroend

## $9 is the complete folder name
Function InstallFolderInternal
  Pop $9
  !define Index 'Line${__LINE__}'
  FindFirst $0 $1 "$9\*"
  StrCmp $0 "" "${Index}-End"
  "${Index}-Loop:"
    StrCmp $1 "" "${Index}-End"
    StrCmp $1 "." "${Index}-Next"
    StrCmp $1 ".." "${Index}-Next"
    IfFileExists "$9\$1\*" 0 "${Index}-Write"
      Push $0
      Push $9
      Push "$9\$1"
      Call InstallFolderInternal
      Pop $9
      Pop $0
      Goto "${Index}-Next"
    "${Index}-Write:"
    FileWrite $UninstallLog "$9\$1$\r$\n"
    "${Index}-Next:"
    FindNext $0 $1
    Goto "${Index}-Loop"
  "${Index}-End:"
  !undef Index
FunctionEnd


; TrimNewlines (copied from NSIS documentation)
; input, top of stack  (e.g. whatever$\r$\n)
; output, top of stack (replaces, with e.g. whatever)
; modifies no other variables.

Function un.TrimNewlines
 Exch $R0
 Push $R1
 Push $R2
 StrCpy $R1 0

 loop:
   IntOp $R1 $R1 - 1
   StrCpy $R2 $R0 1 $R1
   StrCmp $R2 "$\r" loop
   StrCmp $R2 "$\n" loop
   IntOp $R1 $R1 + 1
   IntCmp $R1 0 no_trim_needed
   StrCpy $R0 $R0 $R1

 no_trim_needed:
   Pop $R2
   Pop $R1
   Exch $R0
FunctionEnd

; Function un.RemoveEmptyDirs
  ; Pop $9
  ; !define Index 'Line${__LINE__}'
  ; FindFirst $0 $1 "$INSTDIR$9*"
  ; StrCmp $0 "" "${Index}-End"
  ; "${Index}-Loop:"
    ; StrCmp $1 "" "${Index}-End"
    ; StrCmp $1 "." "${Index}-Next"
    ; StrCmp $1 ".." "${Index}-Next"
      ; Push $0
      ; Push $1
      ; Push $9
      ; Push "$9$1\"
      ; Call un.RemoveEmptyDirs
      ; Pop $9
      ; Pop $1
      ; Pop $0
    ; "${Index}-Remove:"
    ; RMDir "$INSTDIR$9$1"
    ; "${Index}-Next:"
    ; FindNext $0 $1
    ; Goto "${Index}-Loop"
  ; "${Index}-End:"
  ; FindClose $0
  ; !undef Index
; FunctionEnd

