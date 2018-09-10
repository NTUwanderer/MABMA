#ifndef MKPATH_H
#define MKPATH_H
/* Function with behaviour like `mkdir -p'  
	from: http://niallohiggins.com/2009/01/08/mkpath-mkdir-p-alike-in-c-for-unix/
*/
#include <cstdlib>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <libgen.h>
#include <string.h>

int MKPATH(const char *s, mode_t mode);

#endif