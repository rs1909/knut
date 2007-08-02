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
    paramview.h constqtgui.h compthread.h screendialog.h )
QT4_ADD_RESOURCES( ICONS mainwindow.qrc )

ADD_EXECUTABLE( pdde-gui main.cpp
    mainwindow.cpp mainwindow.h 
    paramview.cpp paramview.h
    plotdata.cpp plotdata.h 
    plotwindow.cpp plotwindow.h
    constqtgui.cpp constqtgui.h
    compthread.cpp compthread.h
    screendialog.h
    ${MOC_SRCS} ${ICONS} )

TARGET_LINK_LIBRARIES( pdde-gui
	${PDDE_QTMAIN}
	${QT_QTCORE_LIBRARY}
	${QT_QTGUI_LIBRARY}
	${QT_QTXML_LIBRARY}
	${QT_QTSVG_LIBRARY}
	pddelib laarpack umfpack amd cblas atlas ${CMAKE_DL_LIBS} )

QT4_WRAP_CPP( MOC_CONVERT constqtgui.h )
ADD_EXECUTABLE( cfiletoxml cfiletoxml.cpp constqtgui.cpp constqtgui.h ${MOC_CONVERT} )
TARGET_LINK_LIBRARIES( cfiletoxml 
	${QT_QTCORE_LIBRARY}
	${QT_QTXML_LIBRARY}
	pddelib laarpack umfpack amd cblas atlas ${CMAKE_DL_LIBS} )

IF(WIN32)
  SET_TARGET_PROPERTIES( pdde-gui PROPERTIES LINK_FLAGS "-mthreads -Wl,-enable-stdcall-fixup -Wl,-enable-auto-import -Wl,-enable-runtime-pseudo-reloc -Wl,-s -Wl,-subsystem,windows pdde_res.o" )
  ADD_CUSTOM_COMMAND( TARGET pdde-gui PRE_LINK COMMAND windres -i pdde.rc -o pdde_res.o --include-dir=. )
ENDIF(WIN32)

INSTALL( TARGETS pdde-gui cfiletoxml DESTINATION ${PDDE_BIN_DIR} )
