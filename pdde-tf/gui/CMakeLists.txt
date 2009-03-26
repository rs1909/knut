INCLUDE( ${QT_USE_FILE} )

LINK_DIRECTORIES( ${CBLAS_PATH} ${ATLAS_PATH} )

INCLUDE_DIRECTORIES( ${QT_QTXML_INCLUDE_DIR}
	${QT_QTSVG_INCLUDE_DIR}
	${PDDE_SOURCE_DIR}/src
	${PDDE_BINARY_DIR}/src
	${PDDE_SOURCE_DIR}/ATLAS/include
	${PDDE_SOURCE_DIR}/laarpack
	${PDDE_SOURCE_DIR}/UMFPACK/UMFPACK/Include 
	${PDDE_SOURCE_DIR}/UMFPACK/AMD/Include
	${PDDE_SOURCE_DIR}/UMFPACK/UFconfig )

QT4_WRAP_CPP( MOC_SRCS mainwindow.h plotdata.h plotwindow.h 
    paramview.h constqtgui.h compthread.h screendialog.h)
QT4_ADD_RESOURCES( ICONS mainwindow.qrc )

QT4_WRAP_CPP( MOC_CONVERT constqtgui.h )

ADD_EXECUTABLE( cfiletoxml cfiletoxml.cpp constqtgui.cpp constqtgui.h ${MOC_CONVERT} )
TARGET_LINK_LIBRARIES( cfiletoxml 
	${QT_QTCORE_LIBRARY}
	${QT_QTXML_LIBRARY}
	${COREFOUNDATION_LIBRARY}
	pddelib laarpack umfpack amd ${ATLAS_CBLAS} ${ATLAS_ATLAS} ${CMAKE_DL_LIBS} )

IF(APPLE)
  QT4_WRAP_CPP( MOC_SRCS macopenevent.h)
  
  ADD_EXECUTABLE( PDDE-Cont MACOSX_BUNDLE 
    main.cpp
    mainwindow.cpp mainwindow.h 
    paramview.cpp paramview.h
    plotdata.cpp plotdata.h 
    plotwindow.cpp plotwindow.h
    constqtgui.cpp constqtgui.h
    compthread.h
    screendialog.h
    macopenevent.h
    ${MOC_SRCS} ${ICONS} MacOS/pdde.icns
    ${PDDE_SOURCE_DIR}/src/matrix.h ${PDDE_SOURCE_DIR}/src/pddesys.h)

  SET_SOURCE_FILES_PROPERTIES(
    MacOS/pdde.icns
    PROPERTIES MACOSX_PACKAGE_LOCATION Resources)
  SET_SOURCE_FILES_PROPERTIES(
    ${PDDE_SOURCE_DIR}/src/matrix.h ${PDDE_SOURCE_DIR}/src/pddesys.h
    PROPERTIES MACOSX_PACKAGE_LOCATION Headers)

  TARGET_LINK_LIBRARIES( PDDE-Cont
  	${PDDE_QTMAIN}
  	${QT_QTCORE_LIBRARY}
  	${QT_QTGUI_LIBRARY}
  	${QT_QTXML_LIBRARY}
  	${QT_QTSVG_LIBRARY}
  	${COREFOUNDATION_LIBRARY}
  	pddelib laarpack umfpack amd ${ATLAS_CBLAS} ${ATLAS_ATLAS} ${CMAKE_DL_LIBS} )

  SET_TARGET_PROPERTIES( PDDE-Cont
    PROPERTIES
    MACOSX_BUNDLE_INFO_PLIST ${PDDE_SOURCE_DIR}/gui/MacOS/Info.plist
  	MACOSX_BUNDLE_INFO_STRING "A numerical continuation software"
    MACOSX_BUNDLE_ICON_FILE "pdde.icns"
    MACOSX_BUNDLE_GUI_IDENTIFIER "????"
    MACOSX_BUNDLE_LONG_VERSION_STRING "${PACKAGE_VERSION}"
    MACOSX_BUNDLE_BUNDLE_NAME "PDDE-Cont"
    MACOSX_BUNDLE_SHORT_VERSION_STRING "${PACKAGE_VERSION}"
    MACOSX_BUNDLE_BUNDLE_VERSION "1.0"
    MACOSX_BUNDLE_COPYRIGHT "Copyright © 2002-2009, Robert Szalai")

  INSTALL( TARGETS PDDE-Cont DESTINATION "." )

ENDIF(APPLE)

IF(NOT APPLE)
  ADD_EXECUTABLE( pdde-gui main.cpp
    mainwindow.cpp mainwindow.h 
    paramview.cpp paramview.h
    plotdata.cpp plotdata.h 
    plotwindow.cpp plotwindow.h
    constqtgui.cpp constqtgui.h
    compthread.h
    screendialog.h
    ${MOC_SRCS} ${ICONS} )

  TARGET_LINK_LIBRARIES( pdde-gui
  	${PDDE_QTMAIN}
  	${QT_QTCORE_LIBRARY}
  	${QT_QTGUI_LIBRARY}
  	${QT_QTXML_LIBRARY}
  	${QT_QTSVG_LIBRARY}
  	pddelib laarpack umfpack amd ${ATLAS_CBLAS} ${ATLAS_ATLAS} ${CMAKE_DL_LIBS} )

  INSTALL( TARGETS pdde-gui cfiletoxml DESTINATION ${PDDE_BIN_DIR} )

  IF(WIN32)
    SET_TARGET_PROPERTIES( cfiletoxml PROPERTIES LINK_FLAGS "-Wl,-s" )
    SET_TARGET_PROPERTIES( pdde-gui PROPERTIES LINK_FLAGS "-mthreads -Wl,-enable-stdcall-fixup -Wl,-enable-auto-import -Wl,-enable-runtime-pseudo-reloc -Wl,-s -Wl,-subsystem,windows pdde_res.o" )
    ADD_CUSTOM_COMMAND( TARGET pdde-gui PRE_LINK COMMAND windres -i pdde.rc -o pdde_res.o --include-dir=. )
  ENDIF(WIN32)
ENDIF(NOT APPLE)

