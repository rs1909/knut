/*
 * "$Id: config.h.in 310 2007-09-22 21:00:56Z mike $"
 *
 * Configuration file for Mini-XML, a small XML-like file parsing library.
 *
 * Copyright 2003-2007 by Michael Sweet.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */

/*
 * Include necessary headers...
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <ctype.h>

/*
 * Version number...
 */

#cmakedefine MXML_VERSION	"@MXML_VERSION@"

/*
 * Do we have the snprintf() and vsnprintf() functions?
 */

#cmakedefine HAVE_SNPRINTF
#cmakedefine HAVE_VSNPRINTF

/*
 * Do we have the strXXX() functions?
 */

#cmakedefine HAVE_STRDUP

/*
 * Do we have threading support?
 */

#cmakedefine HAVE_PTHREAD_H

/*
 * Define prototypes for string functions as needed...
 */

#  ifndef HAVE_STRDUP
extern char	*_mxml_strdup(const char *);
#    define strdup _mxml_strdup
#  endif /* !HAVE_STRDUP */

extern char	*_mxml_strdupf(const char *, ...);
extern char	*_mxml_vstrdupf(const char *, va_list);

#  ifndef HAVE_SNPRINTF
extern int	_mxml_snprintf(char *, size_t, const char *, ...);
#    define snprintf _mxml_snprintf
#  endif /* !HAVE_SNPRINTF */

#  ifndef HAVE_VSNPRINTF
extern int	_mxml_vsnprintf(char *, size_t, const char *, va_list);
#    define vsnprintf _mxml_vsnprintf
#  endif /* !HAVE_VSNPRINTF */

/*
 * End of "$Id: config.h.in 310 2007-09-22 21:00:56Z mike $".
 */
