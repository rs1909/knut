INCLUDE( ${QT_USE_FILE} )

LINK_DIRECTORIES( ${CBLAS_PATH} ${ATLAS_PATH} )
SET(INTERNAL_LIBS knutlib laarpack umfpack amd ${VFGEN_LIBS} mxml)
SET(EXTERNAL_LIBS ${GINAC_LDFLAGS} ${BLAS_LIBS}
        ${CMAKE_DL_LIBS} ${KNUT_FORTRAN_RUNTIME})
SET(QT_LIBS ${KNUT_QTMAIN}
        ${QT_QTCORE_LIBRARY}
        ${QT_QTGUI_LIBRARY}
        ${QT_QTSVG_LIBRARY}
        ${COREFOUNDATION_LIBRARY})

INCLUDE_DIRECTORIES( ${QT_QTXML_INCLUDE_DIR}
	${QT_QTSVG_INCLUDE_DIR}
	${KNUT_SOURCE_DIR}/src
	${KNUT_SOURCE_DIR}/vfgen
        ${GINAC_INCLUDE_DIRS}
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
  
  ADD_EXECUTABLE( ${GUI_NAME} MACOSX_BUNDLE 
    main.cpp
    mainwindow.cpp mainwindow.h 
    paramview.cpp paramview.h
    plotdata.cpp plotdata.h 
    plotwindow.cpp plotwindow.h
    constqtgui.cpp constqtgui.h
    compthread.h compthread.cpp
    screendialog.h
    macopenevent.h
    ${MOC_SRCS} ${ICONS} MacOS/knut.icns
    ${KNUT_SOURCE_DIR}/src/matrix.h ${KNUT_SOURCE_DIR}/src/knutsys.h
    ${KNUT_SOURCE_DIR}/src/knerror.h)

  SET_SOURCE_FILES_PROPERTIES(
    MacOS/knut.icns
    PROPERTIES MACOSX_PACKAGE_LOCATION Resources)
  SET_SOURCE_FILES_PROPERTIES(
    ${KNUT_SOURCE_DIR}/src/matrix.h ${KNUT_SOURCE_DIR}/src/knutsys.h
    ${KNUT_SOURCE_DIR}/src/knerror.h
    PROPERTIES MACOSX_PACKAGE_LOCATION Headers)

  TARGET_LINK_LIBRARIES( ${GUI_NAME} ${QT_LIBS} 
	${INTERNAL_LIBS} ${EXTERNAL_LIBS}
	${COREFOUNDATION_LIBRARY})

  SET_TARGET_PROPERTIES( ${GUI_NAME}
    PROPERTIES
    MACOSX_BUNDLE_INFO_PLIST ${KNUT_BINARY_DIR}/gui/MacOS/Info.plist )

  ADD_CUSTOM_COMMAND(TARGET ${GUI_NAME}
    POST_BUILD
    COMMAND
    mkdir -p "${KNUT_BINARY_DIR}/gui/${GUI_NAME}.app/Contents/Frameworks"\;
    mkdir -p "${KNUT_BINARY_DIR}/gui/${GUI_NAME}.app/Contents/Resources/"\;
    cp -HR ${QT_QTGUI_LIBRARY_RELEASE}/Versions/4/Resources/qt_menu.nib "${KNUT_BINARY_DIR}/gui/${GUI_NAME}.app/Contents/Resources/"\;
# Apparently no plugins are required 
#    cp -HR ${QT_PLUGINS_DIR}/imageformats "${KNUT_BINARY_DIR}/gui/${GUI_NAME}.app/Contents/Plugins"\;
    # SETTING UP THE PLUGIN LOCATIONS
    rm -f "${KNUT_BINARY_DIR}/gui/${GUI_NAME}.app/Contents/Resources/qt.conf"\;
    echo "[Paths]" >"${KNUT_BINARY_DIR}/gui/${GUI_NAME}.app/Contents/Resources/qt.conf"\;
    echo "Plugins = Plugins" >"${KNUT_BINARY_DIR}/gui/${GUI_NAME}.app/Contents/Resources/qt.conf"\;
    cmake -P ${KNUT_SOURCE_DIR}/createAppBundle.cmake ${KNUT_BINARY_DIR}\;
)

  INSTALL(DIRECTORY ${KNUT_BINARY_DIR}/gui/${GUI_NAME}.app DESTINATION ${CMAKE_INSTALL_PREFIX}
    USE_SOURCE_PERMISSIONS)
ENDIF(APPLE)

IF(NOT APPLE)
  ADD_EXECUTABLE( ${GUI_NAME} main.cpp
    mainwindow.cpp mainwindow.h 
    paramview.cpp paramview.h
    plotdata.cpp plotdata.h 
    plotwindow.cpp plotwindow.h
    constqtgui.cpp constqtgui.h
    compthread.h compthread.cpp
    screendialog.h
    ${MOC_SRCS} ${ICONS} )

  IF(WIN32)
    TARGET_LINK_LIBRARIES( ${GUI_NAME} ${QT_LIBS}
      ${INTERNAL_LIBS} AUX_BLAS ${EXTERNAL_LIBS})
    SET_TARGET_PROPERTIES( ${GUI_NAME} PROPERTIES LINK_FLAGS "-mthreads -Wl,-enable-stdcall-fixup -Wl,-enable-auto-import -Wl,-enable-runtime-pseudo-reloc -Wl,-s -Wl,-subsystem,windows knut_res.o" )
    ADD_CUSTOM_COMMAND( TARGET ${GUI_NAME} PRE_LINK COMMAND ${CMAKE_RC_COMPILER} -i ${KNUT_SOURCE_DIR}/gui/knut.rc -o knut_res.o --include-dir=${KNUT_SOURCE_DIR}/gui )
  ELSE(WIN32)
    TARGET_LINK_LIBRARIES( ${GUI_NAME} ${QT_LIBS}
      ${INTERNAL_LIBS} ${EXTERNAL_LIBS})
  ENDIF(WIN32)

  INSTALL( TARGETS ${GUI_NAME} DESTINATION ${KNUT_BIN_DIR} )
	
ENDIF(NOT APPLE)
