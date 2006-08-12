LINK_DIRECTORIES( ${CBLAS_PATH} ${ATLAS_PATH} )

FIND_PACKAGE( Qt4 REQUIRED )
INCLUDE( ${QT_USE_FILE} )

INCLUDE_DIRECTORIES( ${QT_QTXML_INCLUDE_DIR}
	${PDDE_SOURCE_DIR}/src
	${PDDE_SOURCE_DIR}/ATLAS/include
	${PDDE_SOURCE_DIR}/UMFPACK/UMFPACK/Include 
	${PDDE_SOURCE_DIR}/UMFPACK/AMD/Include )

QT4_WRAP_CPP( MOC_SRCS mainwindow.h plotdata.h plotwindow.h 
    paramview.h constqtgui.h compthread.h screendialog.h )

ADD_EXECUTABLE( pdde-gui main.cpp
    mainwindow.cpp mainwindow.h 
    paramview.cpp paramview.h
    plotdata.cpp plotdata.h 
    plotwindow.cpp plotwindow.h
    constqtgui.cpp constqtgui.h
    compthread.cpp compthread.h
    screendialog.h
    ${MOC_SRCS} )

TARGET_LINK_LIBRARIES( pdde-gui
	${QT_QTCORE_LIBRARY}
	${QT_QTGUI_LIBRARY}
	${QT_QTXML_LIBRARY}
	pddelib laarpack umfpack amd cblas atlas pddef2c ${CMAKE_DL_LIBS} )

INSTALL_TARGETS( /bin pdde-gui )
