#------------------------------------------------------------------------------
# finding Qt5
#------------------------------------------------------------------------------

# Tell CMake to run moc when necessary:
set(CMAKE_AUTOMOC ON)
# As moc files are generated in the binary dir, tell CMake
# to always look for includes there:
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Widgets finds its own dependencies (QtGui and QtCore).

SET(INTERNAL_LIBS knutlib ${LAARPACK} ${UMFPACK} ${AMD} ${VFGEN_LIBS} ${KNUT_MXML_LIB})
SET(EXTERNAL_LIBS ${BLAS_LIBS} ${CMAKE_DL_LIBS} ${KNUT_FORTRAN_RUNTIME})

INCLUDE_DIRECTORIES( ${KNUT_SOURCE_DIR}/src
  ${KNUT_BINARY_DIR}/src
  ${KNUT_SOURCE_DIR}/laarpack
  ${KNUT_MXML_INC}
  ${UMF_INCLUDE} 
  ${AMD_INCLUDE})

file (GLOB Resources *.qrc)
QT5_ADD_RESOURCES (ResourceSources ${Resources})

IF(APPLE)
  ADD_EXECUTABLE( ${GUI_NAME} MACOSX_BUNDLE 
    main.cpp
    mainwindow.cpp mainwindow.h 
    paramview.cpp paramview.h
#    plotdata.cpp plotdata.h 
#    plotwindow.cpp plotwindow.h
    chartwindow.cpp chartwindow.h
    constqtgui.cpp constqtgui.h
    compthread.h compthread.cpp
    screendialog.h
    datachart.cpp datachart.h
    ${MOC_SRCS} ${ICONS} MacOS/knut.icns
    ${KNUT_SOURCE_DIR}/src/matrix.h ${KNUT_SOURCE_DIR}/src/knutsys.h
    ${KNUT_SOURCE_DIR}/src/knerror.h
    ${ResourceSources})

  SET_SOURCE_FILES_PROPERTIES(
    MacOS/knut.icns
    PROPERTIES MACOSX_PACKAGE_LOCATION Resources)
  SET_SOURCE_FILES_PROPERTIES(
    ${KNUT_SOURCE_DIR}/src/matrix.h ${KNUT_SOURCE_DIR}/src/knutsys.h
    ${KNUT_SOURCE_DIR}/src/knerror.h
    PROPERTIES MACOSX_PACKAGE_LOCATION Headers)

  TARGET_LINK_LIBRARIES( ${GUI_NAME} 
	${INTERNAL_LIBS} ${EXTERNAL_LIBS} ${COREFOUNDATION_LIBRARY} Qt5::Widgets Qt5::PrintSupport Qt5::Svg)

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
#    plotdata.cpp plotdata.h 
#    plotwindow.cpp plotwindow.h
    chartwindow.cpp chartwindow.h
    constqtgui.cpp constqtgui.h
    compthread.h compthread.cpp
    screendialog.h knut.rc
    datachart.cpp datachart.h
    ${MOC_SRCS} ${ICONS} ${ResourceSources})

  IF(WIN32)
    set_target_properties (${GUI_NAME} PROPERTIES WIN32_EXECUTABLE TRUE)
    # need to add 
    TARGET_LINK_LIBRARIES( ${GUI_NAME} ${QT_LIBS}
      ${INTERNAL_LIBS} ${EXTERNAL_LIBS} Qt5::WinMain Qt5::Widgets Qt5::PrintSupport Qt5::Svg Qt5::Charts)
  ELSE(WIN32)
    TARGET_LINK_LIBRARIES( ${GUI_NAME} ${QT_LIBS}
      ${INTERNAL_LIBS} ${EXTERNAL_LIBS} Qt5::Widgets Qt5::PrintSupport Qt5::Svg Qt5::Charts)
  ENDIF(WIN32)

  INSTALL( TARGETS ${GUI_NAME} DESTINATION ${KNUT_BIN_DIR} )
	
ENDIF(NOT APPLE)
