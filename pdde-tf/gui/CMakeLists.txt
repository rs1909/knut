INCLUDE( ${QT_USE_FILE} )

LINK_DIRECTORIES( ${CBLAS_PATH} ${ATLAS_PATH} )

INCLUDE_DIRECTORIES( ${QT_QTXML_INCLUDE_DIR}
	${QT_QTSVG_INCLUDE_DIR}
	${KNUT_SOURCE_DIR}/src
	${KNUT_BINARY_DIR}/src
	${KNUT_SOURCE_DIR}/ATLAS/include
	${KNUT_SOURCE_DIR}/laarpack
	${KNUT_SOURCE_DIR}/UMFPACK/UMFPACK/Include 
	${KNUT_SOURCE_DIR}/UMFPACK/AMD/Include
	${KNUT_SOURCE_DIR}/UMFPACK/UFconfig )

QT4_WRAP_CPP( MOC_SRCS mainwindow.h plotdata.h plotwindow.h 
    paramview.h constqtgui.h compthread.h screendialog.h)
QT4_ADD_RESOURCES( ICONS mainwindow.qrc )

QT4_WRAP_CPP( MOC_CONVERT constqtgui.h )

IF(APPLE)
  QT4_WRAP_CPP( MOC_SRCS macopenevent.h)
  
  ADD_EXECUTABLE( Knut MACOSX_BUNDLE 
    main.cpp
    mainwindow.cpp mainwindow.h 
    paramview.cpp paramview.h
    plotdata.cpp plotdata.h 
    plotwindow.cpp plotwindow.h
    constqtgui.h
    compthread.h
    screendialog.h
    macopenevent.h
    ${MOC_SRCS} ${ICONS} MacOS/knut.icns
    ${KNUT_SOURCE_DIR}/src/matrix.h ${KNUT_SOURCE_DIR}/src/knutsys.h)

  SET_SOURCE_FILES_PROPERTIES(
    MacOS/knut.icns
    PROPERTIES MACOSX_PACKAGE_LOCATION Resources)
  SET_SOURCE_FILES_PROPERTIES(
    ${KNUT_SOURCE_DIR}/src/matrix.h ${KNUT_SOURCE_DIR}/src/knutsys.h
    PROPERTIES MACOSX_PACKAGE_LOCATION Headers)

  TARGET_LINK_LIBRARIES( Knut
  	${KNUT_QTMAIN}
  	${QT_QTCORE_LIBRARY}
  	${QT_QTGUI_LIBRARY}
  	${QT_QTSVG_LIBRARY}
  	${COREFOUNDATION_LIBRARY}
  	knutlib laarpack umfpack amd mxml
  	${ATLAS_CBLAS} ${ATLAS_ATLAS} ${CMAKE_DL_LIBS})

  SET_TARGET_PROPERTIES( Knut
    PROPERTIES
    MACOSX_BUNDLE_INFO_PLIST ${KNUT_SOURCE_DIR}/gui/MacOS/Info.plist
  	MACOSX_BUNDLE_INFO_STRING "A numerical continuation software"
    MACOSX_BUNDLE_ICON_FILE "knut.icns"
    MACOSX_BUNDLE_GUI_IDENTIFIER "????"
    MACOSX_BUNDLE_LONG_VERSION_STRING "${PACKAGE_VERSION}"
    MACOSX_BUNDLE_BUNDLE_NAME "Knut"
    MACOSX_BUNDLE_SHORT_VERSION_STRING "${PACKAGE_VERSION}"
    MACOSX_BUNDLE_BUNDLE_VERSION "1.0"
    MACOSX_BUNDLE_COPYRIGHT "Copyright © 2002-2009, Robert Szalai")

  INSTALL( TARGETS Knut DESTINATION "." )

ENDIF(APPLE)

IF(NOT APPLE)
  ADD_EXECUTABLE( knut-gui main.cpp
    mainwindow.cpp mainwindow.h 
    paramview.cpp paramview.h
    plotdata.cpp plotdata.h 
    plotwindow.cpp plotwindow.h
    constqtgui.h
    compthread.h
    screendialog.h
    ${MOC_SRCS} ${ICONS} )

  TARGET_LINK_LIBRARIES( knut-gui
  	${KNUT_QTMAIN}
  	${QT_QTCORE_LIBRARY}
  	${QT_QTGUI_LIBRARY}
  	${QT_QTSVG_LIBRARY}
  	knutlib laarpack umfpack amd mxml
  	${ATLAS_CBLAS} ${ATLAS_ATLAS} ${CMAKE_DL_LIBS})

  INSTALL( TARGETS knut-gui cfiletoxml DESTINATION ${KNUT_BIN_DIR} )

  IF(WIN32)
    SET_TARGET_PROPERTIES( cfiletoxml PROPERTIES LINK_FLAGS "-Wl,-s" )
    SET_TARGET_PROPERTIES( knut-gui PROPERTIES LINK_FLAGS "-mthreads -Wl,-enable-stdcall-fixup -Wl,-enable-auto-import -Wl,-enable-runtime-pseudo-reloc -Wl,-s -Wl,-subsystem,windows knut_res.o" )
    ADD_CUSTOM_COMMAND( TARGET knut-gui PRE_LINK COMMAND windres -i knut.rc -o knut_res.o --include-dir=. )
  ENDIF(WIN32)
ENDIF(NOT APPLE)
