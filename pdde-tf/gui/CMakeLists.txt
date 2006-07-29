# PROJECT( GUI )

LINK_DIRECTORIES( ${CBLAS_PATH} ${ATLAS_PATH} )

FIND_PACKAGE( Qt4 REQUIRED )
INCLUDE( ${QT_USE_FILE} )

INCLUDE_DIRECTORIES( ${QT_QTXML_INCLUDE_DIR}
	${PDDE_SOURCE_DIR}/src
	${PDDE_SOURCE_DIR}/ATLAS/include
	${PDDE_SOURCE_DIR}/UMFPACK/UMFPACK/Include 
	${PDDE_SOURCE_DIR}/UMFPACK/AMD/Include )

QT4_WRAP_CPP( MOC_SRCS mainwindow.h plotdata.h paramview.h constants.h 
    compthread.h )

ADD_EXECUTABLE( gui main.cpp
    mainwindow.cpp mainwindow.h 
    paramview.cpp paramview.h
    plotdata.cpp plotdata.h 
    constants.cpp constants.h
    compthread.cpp compthread.h
    ${MOC_SRCS} )

TARGET_LINK_LIBRARIES( gui
	${QT_QTCORE_LIBRARY}
	${QT_QTGUI_LIBRARY}
	${QT_QTXML_LIBRARY}
	pddelib laarpack umfpack amd cblas atlas pddef2c ${CMAKE_DL_LIBS} )
