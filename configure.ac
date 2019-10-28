AC_INIT(RPEIF, 1.0.1)

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi

CXX=`"${R_HOME}/bin/R" CMD config CXX`
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`

AC_PROG_CXX

AC_LANG(C++)

AC_CHECK_HEADERS(stdint.h inttypes.h limits.h climits)

# Check if unsigned long long is supported (its not standard yet)
AC_COMPILE_IFELSE([
	AC_LANG_SOURCE([[int main() {  unsigned long long x = 100LL; return 0; }]])
	],
	[
		AC_MSG_RESULT([Type unsigned long long is supported])
		AC_DEFINE([HAVE_UNSIGNED_LONG_LONG], [1], [Define to 1 if unsigned long long is supported])
	],
	[AC_MSG_RESULT([Type unsigned long long is not supported])]
)

# Check if support for uint16_t and uint8_t is complete
if test "$ac_cv_header_stdint_h" == yes; then
	AC_COMPILE_IFELSE([
		AC_LANG_PROGRAM([[#define __STDC_LIMIT_MACROS 1
		#include <stdint.h>]],
					[[uint16_t x = UINT16_MAX; uint8_t y = UINT8_MAX;]])
	],
	[
		AC_MSG_RESULT([uint16_t and unit8_t are available in stdint.h])
		AC_DEFINE([HAVE_UINT8_16_MAX], [1], [Define to 1 if types uint16_t and uint8_t are in stdint.h])
		suitable_integer=yes
	],
	[AC_MSG_RESULT([uint16_t and uint8_t are not available in stdint.h])]
	)
fi
if test "$ac_cv_header_inttypes_h" == yes; then
	AC_COMPILE_IFELSE([
		AC_LANG_PROGRAM([[#define __STDC_LIMIT_MACROS 1
		#include <inttypes.h>]],
					[[uint16_t x = UINT16_MAX; uint8_t y = UINT8_MAX;]])
	],
	[
		AC_MSG_RESULT([uint16_t and unit8_t are available in inttypes.h])
		AC_DEFINE([HAVE_UINT8_16_MAX], [1], [Define to 1 if types uint16_t and uint8_t are in inttypes.h])
		suitable_integer=yes
	],
	[AC_MSG_RESULT([uint16_t and uint8_t are not available in inttypes.h])]
	)
fi

# if test "$suitable_integer" != "yes"; then
# 	AC_MSG_FAILURE([no suitable integral type available])
# fi;

# Check for builtin function to "count trailing zeros" (ctz)
AC_COMPILE_IFELSE([
	AC_LANG_SOURCE([[int main() {  int x = __builtin_ctzl(2UL); return 0; }]])
	],
	[
		AC_MSG_RESULT([Found builtin CTZ method (for unsigned long)])
		AC_DEFINE([HAVE_GCC_CTZL], [1], [Define to 1 if you have function __builtin_ctzl(unsigned long)])
	],
	[AC_MSG_RESULT([No builtin CTZ method (for unsigned long)])]
)

AC_COMPILE_IFELSE([
	AC_LANG_SOURCE([[int main() {  int x = __builtin_ctzll(2ULL); return 0; }]])
	],
	[
		AC_MSG_RESULT([Found builtin CTZ method (for unsigned long long)])
		AC_DEFINE([HAVE_GCC_CTZLL], [1], [Define to 1 if you have function __builtin_ctzll(unsigned long long)])
	],
	[AC_MSG_RESULT([No builtin CTZ method (for unsigned long long)])]
)

AC_CHECK_FUNCS(ffsl)
AC_CHECK_FUNCS(ffsll)

# Check for builtin function to "count set bits" (popcountl)
AC_COMPILE_IFELSE([
	AC_LANG_SOURCE([[int main() {  int x = __builtin_popcountl(10UL); return 0; }]])
],
	[
		AC_MSG_RESULT([Found builtin popcount method (for unsigned long)])
		AC_DEFINE([HAVE_BUILTIN_POPCOUNTL], [1], [Define to 1 if you have function __builtin_popcountl(unsigned long)])
	],
	[AC_MSG_RESULT([No builtin popcount method (for unsigned long)])]
)

AC_COMPILE_IFELSE([
	AC_LANG_SOURCE([[int main() {  int x = __builtin_popcountll(10ULL); return 0; }]])
],
	[
		AC_MSG_RESULT([Found builtin popcount method (for unsigned long long)])
		AC_DEFINE([HAVE_BUILTIN_POPCOUNTLL], [1], [Define to 1 if you have function __builtin_popcountl(unsigned long long)])
	],
	[AC_MSG_RESULT([No builtin popcount method (for unsigned long long)])]
)

# Check for pthreads support
oldCFLAGS="$CFLAGS"
CFLAGS="$CFLAGS -pthread"
AC_CHECK_HEADERS(pthread.h)

CFLAGS="$oldCFLAGS"

AC_SUBST(CXX11FLAGS)

AC_CONFIG_HEADERS([src/autoconfig.h])
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

##
## Check for OpenMP support
##
AC_OPENMP
# since some systems have broken OMP libraries
# we also check that the actual package will work
ac_pkg_openmp_cxx=no
if test -n "${OPENMP_CXXFLAGS}"; then
  AC_MSG_CHECKING([if OpenMP works in C++])
  AC_LANG_CONFTEST(
  [AC_LANG_PROGRAM([[#include <omp.h>]], [[ return omp_get_num_threads (); ]])])
  PKG_CFLAGS="${OPENMP_CXXFLAGS}" PKG_LIBS="${OPENMP_CXXFLAGS}" "$RBIN" CMD SHLIB conftest.c 1>&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD && "$RBIN" --vanilla -q -e "dyn.load(paste('conftest',.Platform\$dynlib.ext,sep=''))" 1>&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD && ac_pkg_openmp_cxx=yes
  AC_MSG_RESULT([${ac_pkg_openmp_cxx}])
fi
# if ${ac_pkg_openmp_cxx} = "yes" then we have OMP, otherwise it will be "no"
if test "${ac_pkg_openmp_cxx}" = no; then
  OPENMP_CXXFLAGS=''
else
    AC_DEFINE([HAVE_OPENMP_C], [1], [Define to 1 if you have a C++ compiler with OpenMP support])
fi
AC_SUBST(OPENMP_CXXFLAGS)