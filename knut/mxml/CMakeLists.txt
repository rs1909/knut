SET( MXML_VERSION "2.5" )
INCLUDE (CheckIncludeFiles)
INCLUDE (CheckFunctionExists)

CHECK_FUNCTION_EXISTS(snprintf HAVE_SNPRINTF)
CHECK_FUNCTION_EXISTS(vsnprintf HAVE_VSNPRINTF)
CHECK_FUNCTION_EXISTS(strdup HAVE_STRDUP)
# there's a memory leak in Mini-XML if we enable pthread support
# CHECK_INCLUDE_FILES(pthread.h HAVE_PTHREAD_H)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config.h)

INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_BINARY_DIR} )

ADD_LIBRARY( mxml 
  mxml-entity.c
  mxml-index.c 
  mxml-private.c mxml-private.h
  mxml-get.c
  mxml-set.c
  mxml-attr.c
  mxml-file.c
  mxml-node.c	
  mxml-search.c	
  mxml-string.c
  mxml.h
  ${CMAKE_CURRENT_BINARY_DIR}/config.h )

SET_TARGET_PROPERTIES( mxml PROPERTIES COMPILE_FLAGS "-D_THREAD_SAFE -D_REENTRANT" )
