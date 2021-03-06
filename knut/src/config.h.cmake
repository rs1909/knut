/* Define to the full name of this package. */
#define PACKAGE_NAME "@PACKAGE_NAME@"

/* Define to the version of this package. */
#define PACKAGE_VERSION "@PACKAGE_VERSION@"

/* "SVN REVISION" */
#define PACKAGE_REVISION "@PACKAGE_REVISION@"

#define PACKAGE_COPYRIGHT "@PACKAGE_COPYRIGHT@"
#define PACKAGE_URL "@PACKAGE_URL@"

/* definitions for compiling a system definition file */
#define CMAKE_CXX_COMPILER                  "@CMAKE_CXX_COMPILER@"
#define CMAKE_CXX_FLAGS                     "@CMAKE_CXX_FLAGS@"
#define CMAKE_SHARED_LIBRARY_C_FLAGS        "@CMAKE_SHARED_LIBRARY_C_FLAGS@"
#define CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "@CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS@"
#define KNUT_NO_EXCEPTIONS                  "@KNUT_NO_EXCEPTIONS@"
#define CMAKE_INSTALL_PREFIX                "@CMAKE_INSTALL_PREFIX@"
#define CMAKE_OSX_ARCHITECTURES             "@CMAKE_OSX_ARCHITECTURES@"
#define KNUT_INCLUDE_DIR                    "@KNUT_INCLUDE_DIR@"
#define KNUT_SOURCE_DIR                     "@KNUT_SOURCE_DIR@"

/* whether build the GUI */
#cmakedefine KNUT_GUI
#cmakedefine GINAC_FOUND

/* whether 64bit BLAS installed */
#cmakedefine LONGBLAS
