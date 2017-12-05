#------------------------------------------------------------------------------
# finding Qt5
#------------------------------------------------------------------------------

# Tell CMake to run moc when necessary:
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
# As moc files are generated in the binary dir, tell CMake
# to always look for includes there:
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Widgets finds its own dependencies (QtGui and QtCore).

SET(INTERNAL_LIBS knutlib ${LAARPACK} ${UMFPACK} ${AMD} ${VFGEN_LIBS} ${KNUT_MXML_LIB})
SET(EXTERNAL_LIBS ${BLAS_LIBS} ${CMAKE_DL_LIBS} ${KNUT_FORTRAN_RUNTIME})
set(GUI_SOURCE_FILES main.cpp
    mainwindow.cpp mainwindow.h 
    paramview.cpp paramview.h
    plotdata.cpp plotdata.h 
    plotwindow.cpp plotwindow.h
#    chartwindow.cpp chartwindow.h
    constqtgui.cpp constqtgui.h
    compthread.h compthread.cpp
    screendialog.h
#    datachart.cpp datachart.h 
    mainwindow.qrc
)

INCLUDE_DIRECTORIES( ${KNUT_SOURCE_DIR}/src
  ${KNUT_BINARY_DIR}/src
  ${KNUT_SOURCE_DIR}/laarpack
  ${KNUT_MXML_INC}
  ${UMF_INCLUDE} 
  ${AMD_INCLUDE})

IF(APPLE)
  ADD_EXECUTABLE( ${GUI_NAME} MACOSX_BUNDLE 
    ${GUI_SOURCE_FILES}
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

  IF(WIN32)
    ADD_EXECUTABLE( ${GUI_NAME} WIN32
        ${GUI_SOURCE_FILES}
        ${MOC_SRCS} ${ICONS} ${ResourceSources})
    set_target_properties (${GUI_NAME} PROPERTIES WIN32_EXECUTABLE TRUE)
    # need to add 
    TARGET_LINK_LIBRARIES( ${GUI_NAME} ${QT_LIBS}
      ${INTERNAL_LIBS} ${EXTERNAL_LIBS} Qt5::WinMain Qt5::Widgets Qt5::PrintSupport Qt5::Svg)
  ELSE()
    ADD_EXECUTABLE( ${GUI_NAME}
        ${GUI_SOURCE_FILES}
        ${MOC_SRCS} ${ICONS} ${ResourceSources})
    TARGET_LINK_LIBRARIES( ${GUI_NAME} ${QT_LIBS}
      ${INTERNAL_LIBS} ${EXTERNAL_LIBS} Qt5::Widgets Qt5::PrintSupport Qt5::Svg)
  ENDIF()

  INSTALL( TARGETS ${GUI_NAME} DESTINATION ${KNUT_BIN_DIR} )
	
ENDIF(NOT APPLE)

set_target_properties(${GUI_NAME} PROPERTIES CXX_STANDARD 17)
