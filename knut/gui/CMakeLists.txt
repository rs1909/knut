#------------------------------------------------------------------------------
# finding Qt5
#------------------------------------------------------------------------------

# Tell CMake to run moc when necessary:
set(CMAKE_AUTOMOC ON)
# As moc files are generated in the binary dir, tell CMake
# to always look for includes there:
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Widgets finds its own dependencies (QtGui and QtCore).

SET(INTERNAL_LIBS knutlib laarpack umfpack amd ${VFGEN_LIBS} mxml)
SET(EXTERNAL_LIBS ${GINAC_LDFLAGS} ${BLAS_LIBS}
        ${CMAKE_DL_LIBS} ${KNUT_FORTRAN_RUNTIME} ${COREFOUNDATION_LIBRARY})

INCLUDE_DIRECTORIES( ${KNUT_SOURCE_DIR}/src
	${KNUT_SOURCE_DIR}/vfgen
        ${GINAC_INCLUDE_DIRS}
	${KNUT_BINARY_DIR}/src
	${KNUT_SOURCE_DIR}/ATLAS/include
	${KNUT_SOURCE_DIR}/laarpack
	${KNUT_SOURCE_DIR}/UMFPACK/UMFPACK/Include 
	${KNUT_SOURCE_DIR}/UMFPACK/AMD/Include
	${KNUT_SOURCE_DIR}/UMFPACK/UFconfig )

file (GLOB Resources *.qrc)
QT5_ADD_RESOURCES (ResourceSources ${Resources}) 

IF(APPLE)
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
    ${KNUT_SOURCE_DIR}/src/knerror.h
    ${ResourceSources})

    qt5_use_modules(${GUI_NAME} Widgets)
    qt5_use_modules(${GUI_NAME} PrintSupport)
    qt5_use_modules(${GUI_NAME} Svg)

  SET_SOURCE_FILES_PROPERTIES(
    MacOS/knut.icns
    PROPERTIES MACOSX_PACKAGE_LOCATION Resources)
  SET_SOURCE_FILES_PROPERTIES(
    ${KNUT_SOURCE_DIR}/src/matrix.h ${KNUT_SOURCE_DIR}/src/knutsys.h
    ${KNUT_SOURCE_DIR}/src/knerror.h
    PROPERTIES MACOSX_PACKAGE_LOCATION Headers)

  TARGET_LINK_LIBRARIES( ${GUI_NAME} 
	${INTERNAL_LIBS} ${EXTERNAL_LIBS})

  SET_TARGET_PROPERTIES( ${GUI_NAME}
    PROPERTIES
    MACOSX_BUNDLE_INFO_PLIST ${KNUT_BINARY_DIR}/gui/MacOS/Info.plist )
  ADD_CUSTOM_COMMAND(TARGET ${GUI_NAME}
    POST_BUILD
    COMMAND
    macdeployqt "${KNUT_BINARY_DIR}/gui/${GUI_NAME}.app" )

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
    screendialog.h knut.rc
    ${MOC_SRCS} ${ICONS} ${ResourceSources})

    qt5_use_modules(${GUI_NAME} Widgets)
    qt5_use_modules(${GUI_NAME} PrintSupport)
    qt5_use_modules(${GUI_NAME} Svg)

  IF(WIN32)
    SET(CMAKE_RC_COMPILER "windres")
    SET(CMAKE_RC_COMPILE_OBJECT "<CMAKE_RC_COMPILER> <FLAGS> <DEFINES> -O coff -o <OBJECT> <SOURCE>")
    ENABLE_LANGUAGE(RC) 
    # need to add 
    TARGET_LINK_LIBRARIES( ${GUI_NAME} ${QT_LIBS}
    ${INTERNAL_LIBS} ${EXTERNAL_LIBS} Qt5::WinMain )
#     SET_TARGET_PROPERTIES( ${GUI_NAME} PROPERTIES LINK_FLAGS "-mthreads -Wl,-enable-stdcall-fixup -Wl,-enable-auto-import -Wl,-enable-runtime-pseudo-reloc -Wl,-s -Wl,-subsystem,windows knut_res.o" )
#     ADD_CUSTOM_COMMAND( TARGET ${GUI_NAME} PRE_LINK COMMAND ${CMAKE_RC_COMPILER} -i ${KNUT_SOURCE_DIR}/gui/knut.rc -o knut_res.o --include-dir=${KNUT_SOURCE_DIR}/gui )
  ELSE(WIN32)
    TARGET_LINK_LIBRARIES( ${GUI_NAME} ${QT_LIBS}
      ${INTERNAL_LIBS} ${EXTERNAL_LIBS} )
  ENDIF(WIN32)

  INSTALL( TARGETS ${GUI_NAME} DESTINATION ${KNUT_BIN_DIR} )
	
ENDIF(NOT APPLE)
