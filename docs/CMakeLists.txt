

set(SRCDIR ${CMAKE_CURRENT_SOURCE_DIR})
set(BLDDIR ${CMAKE_CURRENT_BINARY_DIR})

# --------------------------
#  Need sphinx
# --------------------------
find_program(Sphinx_EXECUTABLE
  sphinx-build
)
if (Sphinx_EXECUTABLE)
  set(SPHINX_FOUND 1 CACHE BOOL "Found Sphinx binary")
  message(STATUS "Sphinx_EXECUTABLE found.")
  message(STATUS "Sphinx_EXECUTABLE  = ${Sphinx_EXECUTABLE}")
else ()
  message(STATUS "Sphinx_EXECUTABLE NOT found.")
endif ()

# Enable setting an option
# The -d flag sets the cache directory to be non-hidden
set(Sphinx_OPTS      "" CACHE STRING "Options to be passed to the sphinx executable")

# Handle the build issues
include(SciSphinxFunctions.cmake)

# Will need to figure this out
set(MATHJAXJS "MathJax/MathJax.js?config=TeX-AMS-MML_SVG-full.js")


#---------------------
# Create documentation 
#---------------------
set(RST_FILE_BASE index)
file(GLOB_RECURSE EFITDOC_DEPS RELATIVE ${SRCDIR} *.rst */*.rst)
set(EFITDOC_DEPS ${EFITDOC_DEPS} ../INSTALL.rst ../LICENSE.rst)

# Copy stuff over from docs because sphinxbuild doesn't allow full paths
configure_file(${SRCDIR}/conf.py.in ${BLDDIR}/conf.py)  # 1 copy
configure_file(${SRCDIR}/build_index.html.in ${BLDDIR}/build_index.html)  # 1 copy

# Pull together all ReST files for reference manual
file(GLOB_RECURSE EFITDOC_DEPS RELATIVE ${SRCDIR} *.rst */*.rst */*/*.rst */*/*/*.rst)

# This sets the targets for each build (html, pdf, ...)
  SciSphinxTarget(TARGET docs
                  RST_FILE_BASE ${RST_FILE_BASE}
                  FILE_DEPS ${EFITDOC_DEPS}
                  INSTALL_SUPERDIR doc
                  SPHINX_DOCTREE_DIR "-d${BLDDIR}/doctrees"
                  SPHINX_ADDL_OPTS ${SPHINX_MATHARG}
  )


  # Add top level targets
add_custom_target(docs ALL DEPENDS docs-html)
