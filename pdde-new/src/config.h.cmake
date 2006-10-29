/* Define to the full name of this package. */
#define PACKAGE_NAME "@PACKAGE_NAME@"

/* Define to the version of this package. */
#define PACKAGE_VERSION "@PACKAGE_VERSION@"

/* "SVN REVISION" */
#define PACKAGE_REVISION "@PACKAGE_REVISION@"

/* definitions for compiling a system definition file */
#define CMAKE_CXX_COMPILER                  "@CMAKE_CXX_COMPILER@"
#define CMAKE_CXX_FLAGS                     "@CMAKE_CXX_FLAGS_RELEASE@"
#define CMAKE_SHARED_LIBRARY_C_FLAGS        "@CMAKE_SHARED_LIBRARY_C_FLAGS@"
#define CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "@CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS@"
#define CMAKE_INSTALL_PREFIX                "@CMAKE_INSTALL_PREFIX@"

/* whether build the GUI */
#cmakedefine PDDE_GUI
