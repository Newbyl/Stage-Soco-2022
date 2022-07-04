//******************************************************************************
//	Framework.h
//
//	To include header files based on the targeted operating system.
//******************************************************************************
#ifndef FRAMEWORK_H_
#define FRAMEWORK_H_

#ifdef _WINDOWS
#define WIN32_LEAN_AND_MEAN		// Exclude unusual headers
#include <stdio.h>
#include <stdlib.h>
#include <windows.h>
#include <stdint.h>
#define EXPORT extern "C" __declspec(dllexport)
#define IMPORT extern "C" __declspec(dllimport)
#define PRIVATE
#define TRACE(x) OutputDebugStringA(x);
#endif /* WIN32 */

#ifdef _LINUX
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <dirent.h>
#include <dlfcn.h>
#include <stdint.h>
#include <unistd.h>
#include <link.h>
#define __max(a,b) (((a) > (b)) ? (a) : (b))
#define __min(a,b) (((a) < (b)) ? (a) : (b))
typedef void* HMODULE;
#define WINAPI
#define EXPORT extern "C" __attribute__((visibility("default")))
#define IMPORT
#define PRIVATE extern "C" __attribute__((visibility("hidden")))
#define TRACE(x) printf(x)
#define strcmpi(x,y) strcasecmp((x),(y))
#define strnicmp(x,y,z) strncasecmp((x),(y),(z))
#endif /* _LINUX */

#ifdef __APPLE__
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <dirent.h>
#include <dlfcn.h>
#include <stdint.h>
#include <unistd.h>
#define __max(a,b) (((a) > (b)) ? (a) : (b))
#define __min(a,b) (((a) < (b)) ? (a) : (b))
typedef void* HMODULE;
#define WINAPI
#define EXPORT extern "C" __attribute__((visibility("default")))
#define IMPORT
#define PRIVATE extern "C" __attribute__((visibility("hidden")))
#define TRACE(x) printf(x)
#define strcmpi(x,y) strcasecmp((x),(y))
#define _strcmpi(x,y) strcasecmp((x),(y))
#define strnicmp(x,y,z) strncasecmp((x),(y),(z))
#define _strnicmp(x,y,z) strncasecmp((x),(y),(z))
#endif /* _LINUX */

#endif /* FRAMEWORK_H_ */
