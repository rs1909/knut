;--------------------------------
;Include Modern UI

  !include "MUI.nsh"

;--------------------------------
;General

  !define PACKAGE_NAME "@PACKAGE_NAME@"
  !define PACKAGE_VERSION "@PACKAGE_VERSION@"
  ;Name and file
  Name "${PACKAGE_NAME}"
  OutFile "${PACKAGE_NAME}-${PACKAGE_VERSION}.exe"

  ;Default installation folder
  InstallDir "$PROGRAMFILES\${PACKAGE_NAME}"
  
  ;Get installation folder from registry if available
  InstallDirRegKey HKCU "Software\${PACKAGE_NAME}" ""

;--------------------------------
;Variables

  Var MUI_TEMP
  Var STARTMENU_FOLDER

;--------------------------------
;Interface Settings

  !define MUI_ABORTWARNING

;--------------------------------
;Pages

  !insertmacro MUI_PAGE_LICENSE "COPYING"
  !insertmacro MUI_PAGE_DIRECTORY
  
  ;Start Menu Folder Page Configuration
  !define MUI_STARTMENUPAGE_REGISTRY_ROOT "HKCU" 
  !define MUI_STARTMENUPAGE_REGISTRY_KEY "Software\${PACKAGE_NAME}" 
  !define MUI_STARTMENUPAGE_REGISTRY_VALUENAME "Start Menu Folder"
  
  !insertmacro MUI_PAGE_STARTMENU Application $STARTMENU_FOLDER
  
  !insertmacro MUI_PAGE_INSTFILES
  
  !insertmacro MUI_UNPAGE_CONFIRM
  !insertmacro MUI_UNPAGE_INSTFILES

;--------------------------------
;Languages
 
  !insertmacro MUI_LANGUAGE "English"

;--------------------------------
;Installer Sections

Section "Dummy Section" SecDummy

  SetOutPath "$INSTDIR"
  
  ;ADD YOUR OWN FILES HERE...
  File "COPYING"
  File /r "@CMAKE_INSTALL_PREFIX@\bin"
  File /r "@CMAKE_INSTALL_PREFIX@\include"
  File /r "@CMAKE_INSTALL_PREFIX@\demo"
  File /r "@CMAKE_INSTALL_PREFIX@\matlab"
  
  ;Store installation folder
  WriteRegStr HKCU "Software\${PACKAGE_NAME}" "" $INSTDIR
  
  ;Create uninstaller
  WriteUninstaller "$INSTDIR\Uninstall.exe"
  
  !insertmacro MUI_STARTMENU_WRITE_BEGIN Application
    
    ;Create shortcuts
    CreateDirectory "$SMPROGRAMS\$STARTMENU_FOLDER"
    CreateShortCut "$SMPROGRAMS\$STARTMENU_FOLDER\pdde-gui.lnk" "$INSTDIR\bin\pdde-gui.exe"
    CreateShortCut "$SMPROGRAMS\$STARTMENU_FOLDER\Uninstall.lnk" "$INSTDIR\Uninstall.exe"
  
  !insertmacro MUI_STARTMENU_WRITE_END

SectionEnd

;--------------------------------
;Descriptions
;
;  ;Language strings
;  LangString DESC_SecDummy ${LANG_ENGLISH} "A test section."
;
;  ;Assign language strings to sections
;  !insertmacro MUI_FUNCTION_DESCRIPTION_BEGIN
;    !insertmacro MUI_DESCRIPTION_TEXT ${SecDummy} $(DESC_SecDummy)
;  !insertmacro MUI_FUNCTION_DESCRIPTION_END
; 
;--------------------------------
;Uninstaller Section

Section "Uninstall"

  ;ADD YOUR OWN FILES HERE...
  RMdir /r "$INSTDIR\bin"
  RMdir /r "$INSTDIR\include"
  RMdir /r "$INSTDIR\matlab"
  RMdir /r "$INSTDIR\demo"
  Delete "$INSTDIR\COPYING"

  Delete "$INSTDIR\Uninstall.exe"

  RMDir "$INSTDIR"
  
  !insertmacro MUI_STARTMENU_GETFOLDER Application $MUI_TEMP
    
  Delete "$SMPROGRAMS\$MUI_TEMP\Uninstall.lnk"
  Delete "$SMPROGRAMS\$MUI_TEMP\pdde-gui.lnk"
  
  ;Delete empty start menu parent diretories
  StrCpy $MUI_TEMP "$SMPROGRAMS\$MUI_TEMP"
 
  startMenuDeleteLoop:
	ClearErrors
    RMDir $MUI_TEMP
    GetFullPathName $MUI_TEMP "$MUI_TEMP\.."
    
    IfErrors startMenuDeleteLoopDone
  
    StrCmp $MUI_TEMP $SMPROGRAMS startMenuDeleteLoopDone startMenuDeleteLoop
  startMenuDeleteLoopDone:

  DeleteRegKey /ifempty HKCU "Software\${PACKAGE_NAME}"

SectionEnd