# - SciSphinxFunctions:
# Useful functions for simplifying the setting up of sphinx targets
#
# All functions assume that FindSciSphinx was used and the following are
# defined:
#   Sphinx_EXECUTABLE     = The path to the sphinx command.
#   Sphinx_OPTS           = Options for sphinx
#
# The SciSphinxTarget automates the creation of build and install
# targets.  The make targets are not added to all, and the install's are
# optional.  To make the install options work with install, use the
# add_dependencies command.  For example:
#  add_dependencies(install install-userdocs)
#

#################################################################
#
# SciSphinxFunction
#
# @version $Rev: 1809 $ $Date: 2020-11-21 10:14:57 -0700 (Sat, 21 Nov 2020) $
#
# Copyright &copy; 2012-2020, Tech-X Corporation, Boulder, CO.
# See LICENSE file (EclipseLicense.txt) for conditions of use.
#
#################################################################

include(CMakeParseArguments)

# SciSphinxTarget.cmake
# Automate the defining of the CMake targets
# Args:
#   TARGET:        Target basename.  Targets will be ${TARGET}-<build>,
#                    where <build> is one of html, latex, or pdf.
#   RST_FILE_BASE: Root name of Latex file.  From conf.py
#   SOURCE_DIR:    Directory containing the index.rst.  Defaults
#                    to CMAKE_CURRENT_SOURCE_DIR
#   SPHINX_ADDL_OPTS:   Additional options to Sphinx
#   SPHINX_DOCTREE_ARG: Select cache directory (default .doctrees)
#   FILE_DEPS:      Files that are the dependencies.
#   SPHINX_BUILDS:  Which builds to include.  Default is "html latex pdf"
#                     Possible choices are "html latex pdf singlehtml man"
#   SPHINX_INSTALLS:     Which builds to install.  Default is same as builds
#   NOWARN_NOTMATCH_DIR: Do not warn if file base does not match install dir
#   INSTALL_SUPERDIR:    Name of installation directory up to this one.
#                          Should not be absolute (not include prefix).
#                          Overridden by INSTALL_SUBDIR.
#   INSTALL_SUBDIR:   Name of this subdir for installation.
#                     Should not be absolute (not include prefix).
#
macro(SciSphinxTarget)

# Parse out the args
  set(opts DEBUG;NOWARN_NOTMATCH_DIR) # no-value args
  set(oneValArgs RST_FILE_BASE;TARGET;SPHINX_ADDL_OPTS;SPHINX_DOCTREE_ARG;
    SOURCE_DIR;INSTALL_SUPERDIR;INSTALL_SUBDIR)
  set(multValArgs FILE_DEPS;ALL_BUILDS) # e.g., lists
  cmake_parse_arguments(FD "${opts}" "${oneValArgs}" "${multValArgs}" ${ARGN})
#
# Defaults
#
  if (NOT DEFINED FD_SOURCE_DIR)
    set(FD_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
  endif ()
  if (DEFINED FD_SPHINX_BUILDS)
# If there is a pdf build, there must be a latex build
    list(FIND FD_SPHINX_BUILDS pdf pdfindx)
    if (NOT pdfindx EQUAL -1)
      list(FIND FD_SPHINX_BUILDS latex latexindx)
      if (latexindx EQUAL -1)
        list(APPEND FD_SPHINX_BUILDS latex)
      endif ()
    endif ()
  else ()
    set(FD_SPHINX_BUILDS)
    list(APPEND FD_SPHINX_BUILDS html latex pdf)
  endif ()
  if (NOT DEFINED FD_SPHINX_INSTALLS)
    set(FD_SPHINX_INSTALLS ${FD_SPHINX_BUILDS})
  endif ()
  if (FD_INSTALL_SUBDIR)
    set(instdir ${FD_INSTALL_SUBDIR})
  elseif (FD_INSTALL_SUPERDIR)
    set(instdir ${FD_INSTALL_SUPERDIR}/${thissubdir})
  else ()
    set(instdir ${CMAKE_INSTALL_PREFIX})
  endif ()
  if (NOT DEFINED USE_WHOOSH )
    set(USE_WHOOSH false)
  endif ()
  # if (NOT DEFINED FD_LATEX_NUMTIMES)
    # set(FD_LATEX_NUMTIMES 1)
  # endif ()

#
#  Basic sanity checks
#
  get_filename_component(thissubdir ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  if (NOT NOWARN_NOTMATCH_DIR)
    set(WARN_NOTMATCH_DIR)
  endif ()
  if (WARN_NOTMATCH_DIR)
    if (NOT "${thissubdir}" STREQUAL "${FD_RST_FILE_BASE}")
      message(WARNING "Main rst file base, ${FD_RST_FILE_BASE}, does not match subdirectory name, ${thissubdir}.")
    endif ()
  endif ()
  if (NOT DEFINED FD_TARGET)
    message(WARNING "SciSphinxTarget called without TARGET.")
    return()
  endif ()
  if (NOT DEFINED FD_FILE_DEPS)
    message(WARNING "SciSphinxTarget called without FILE_DEPS.")
    return()
  endif ()
  if (NOT DEFINED FD_RST_FILE_BASE)
    message(WARNING "SciSphinxTarget called without RST_FILE_BASE.")
    return()
  endif ()
  if (NOT DEFINED Sphinx_EXECUTABLE)
    message(WARNING "SciSphinxTarget called without Sphinx_EXECUTABLE.")
    return()
  endif ()
  if (FD_DEBUG)
    message(STATUS "
")
    message("--------- SciSphinxTarget defining targets for ${FD_TARGET} ---------")
    message(STATUS "[SciSphinxFunctions]: TARGET = ${FD_TARGET} ")
    message(STATUS "[SciSphinxFunctions]: RST_FILE_BASE = ${FD_RST_FILE_BASE} ")
    message(STATUS "[SciSphinxFunctions]: Sphinx_EXECUTABLE = ${Sphinx_EXECUTABLE} ")
    message(STATUS "[SciSphinxFunctions]: Sphinx_OPTS = ${Sphinx_OPTS} ")
    message(STATUS "[SciSphinxFunctions]: SPHINX_ADDL_OPTS = ${FD_SPHINX_ADDL_OPTS} ")
    message(STATUS "[SciSphinxFunctions]: SPHINX_DOCTREE_ARG = ${FD_SPHINX_DOCTREE_ARG} ")
  endif ()
  if (FD_SPHINX_DOCTREE_ARG)
    message(WARNING "SPHINX_DOCTREE_ARG set.  Different for each build?")
  endif ()

#
#  Set the output for the standard builds
#
  if (NOT EXISTS ${FD_SOURCE_DIR}/${FD_RST_FILE_BASE}.rst)
    set(html_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/html/${FD_RST_FILE_BASE}.html)
  else ()
    set(html_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/html/index.html)
  endif ()
  set(singlehtml_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/singlehtml/index.html)
  set(latex_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/pdf/${FD_RST_FILE_BASE}.tex)
  set(pdf_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/pdf/${FD_RST_FILE_BASE}.pdf)
  set(man_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/man/index.man)

  foreach (build ${FD_SPHINX_BUILDS})
    set(BUILD_DEPS)
    set(TOUCH_CMD_ARGS)
# touch to redo bibs
    if (SCI_BIBTEX_REFS)
     find_program(TOUCH touch)
      if (TOUCH)
        set(TOUCH_CMD_ARGS COMMAND touch ARGS ${SCI_BIBTEX_REFS})
      else ()
# Could do this better.  See mklinks.sh
        message(WARNING "touch not known on this platform.  bibs might be out of date.")
      endif ()
    endif ()
    set(${build}_DIR ${CMAKE_CURRENT_BINARY_DIR}/${build})
# Latex is actually for pdf which is below
    if (${build} STREQUAL latex)
      set(${build}_DIR ${CMAKE_CURRENT_BINARY_DIR}/pdf)
    endif ()

# There is something weird about passing blank spaces into COMMAND
# so this method fixes the problems that arise if Sphinx_OPTS is not defined
    set(all_opts -b ${build} -c ${CMAKE_CURRENT_BINARY_DIR} ${Sphinx_OPTS} ${FD_SPHINX_ADDL_OPTS} ${FD_SPHINX_DOCTREE_ARG})
    if (NOT ${build} STREQUAL pdf)
      if (${build} STREQUAL html AND MathJax_ROOT_DIR AND NOT "${MATHJAXJS}" MATCHES "^http")
# If html build and have mathjax and not remote mathjax,
# create a target to make the directory down to the mathjax config dir,
# copy in mathjaxconf.js, and then copy in MathJax.
# The output is mathjaxconf.js in the MathJax config dir, so this
# will not be redone if that file is there.
# Also define BUILD_DEPS, a dependency to be added to the build of the
# html docs so that this gets done when that target is invoked.
        set(STATIC_DIR ${CMAKE_CURRENT_BINARY_DIR}/html/_static)
        set(MathJax_DIR ${STATIC_DIR}/MathJax)
        set(MathJax_CONFIG_DIR ${MathJax_DIR}/config)
        set(MathJax_OUTPUT ${MathJax_CONFIG_DIR}/mathjaxconf.js)
        add_custom_command(
            OUTPUT ${MathJax_OUTPUT}
            COMMAND ${CMAKE_COMMAND}
            ARGS -E make_directory ${MathJax_CONFIG_DIR}
            COMMAND ${CMAKE_COMMAND}
            ARGS -E copy ${CMAKE_SOURCE_DIR}/customizations/mathjaxconf.js ${MathJax_CONFIG_DIR}
            COMMAND ${CMAKE_COMMAND}
            ARGS -E copy_directory ${MathJax_ROOT_DIR} ${MathJax_DIR}
            DEPENDS ${CMAKE_SOURCE_DIR}/customizations/mathjaxconf.js
        )
        add_custom_target(${FD_TARGET}-mjconf DEPENDS ${MathJax_OUTPUT})
        set(BUILD_DEPS ${MathJax_OUTPUT})
      endif ()
      add_custom_command(
        OUTPUT ${${build}_OUTPUT}
        COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/gen_doxyrst.py
        ARGS ${CMAKE_SOURCE_DIR}/efit
        COMMAND ${DOXYGEN_EXECUTABLE}
        #COMMAND ln -sf ${CMAKE_CURRENT_BINARY_DIR}/xml ${CMAKE_CURRENT_SOURCE_DIR}
        COMMAND ${Sphinx_EXECUTABLE}
        ARGS ${all_opts} ${FD_SOURCE_DIR} ${${build}_DIR}
        ${TOUCH_CMD_ARGS}
        COMMAND ${Sphinx_EXECUTABLE}
        ARGS ${all_opts} ${FD_SOURCE_DIR} ${${build}_DIR}
        DEPENDS ${FD_FILE_DEPS}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      )
# Add the above found build deps, if any as a dependency to the target.
      # message(STATUS "${FD_TARGET}-${build} DEPENDS ${BUILD_DEPS} ${${build}_OUTPUT}")
      add_custom_target(${FD_TARGET}-${build} DEPENDS ${BUILD_DEPS} ${${build}_OUTPUT})
      if (USE_WHOOSH AND "${build}" STREQUAL html)
        message(STATUS "XXXXXXXXXXX")
        set(whoosh_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/html/whoosh.txt)
        set(all_opts -b whoosh -c ${CMAKE_CURRENT_BINARY_DIR} ${Sphinx_OPTS} ${FD_SPHINX_ADDL_OPTS} ${FD_SPHINX_DOCTREE_ARG})
        add_custom_command(
          OUTPUT ${whoosh_OUTPUT}
          COMMAND ${Sphinx_EXECUTABLE}
          ARGS ${all_opts} ${FD_SOURCE_DIR} ${${build}_DIR}
          DEPENDS ${FD_FILE_DEPS}
        )
        add_custom_target(${FD_TARGET}-whoosh DEPENDS ${whoosh_OUTPUT})
      endif ()
    endif ()
  endforeach ()

#
#  PDF is special
#   This must be make, as sphinx generates a unix makefile
#
  add_custom_command(
    OUTPUT ${pdf_OUTPUT}
    COMMAND make all-pdf LATEXMKOPTS=-bibtex
    DEPENDS ${latex_OUTPUT}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/pdf
  )
  add_custom_target(${FD_TARGET}-pdf DEPENDS ${pdf_OUTPUT})

#
#  Each install is a one-off
#
  list(FIND FD_SPHINX_INSTALLS "pdf" indx)
  if (NOT indx EQUAL -1)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/pdf/${FD_RST_FILE_BASE}.pdf
      DESTINATION "${instdir}"
      PERMISSIONS OWNER_WRITE OWNER_READ GROUP_WRITE GROUP_READ WORLD_READ
      OPTIONAL
    )
  endif ()
  list(FIND FD_SPHINX_INSTALLS "html" indx)
  if (NOT indx EQUAL -1)
    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html
      DESTINATION ${instdir}
      FILE_PERMISSIONS OWNER_WRITE OWNER_READ GROUP_WRITE GROUP_READ WORLD_READ
      OPTIONAL
      PATTERN ".buildinfo" EXCLUDE
      PATTERN "CMakeLists.txt" EXCLUDE
      PATTERN "README" EXCLUDE
# The below not found in */*.html
      PATTERN "ajax-loader.gif" EXCLUDE
      PATTERN "bootstrap-2.3.2" EXCLUDE
      PATTERN "bootstrap-3.3.4" EXCLUDE
      PATTERN "bootswatch-2.3.2" EXCLUDE
      PATTERN "bootswatch-3.3.4" EXCLUDE
      PATTERN "banner_logo_800_text_glow.jpg" EXCLUDE
      PATTERN "comment-bright.png" EXCLUDE
      PATTERN "comment-close.png" EXCLUDE
      PATTERN "comment.png" EXCLUDE
      PATTERN "down-pressed.png" EXCLUDE
      PATTERN "down.png" EXCLUDE
      PATTERN "favicon.ico" EXCLUDE
      PATTERN "file.png" EXCLUDE
      PATTERN "minus.png" EXCLUDE
      PATTERN "plus.png" EXCLUDE
      PATTERN "templates" EXCLUDE
      PATTERN "up-pressed.png" EXCLUDE
      PATTERN "up.png" EXCLUDE
    )
  endif ()
  list(FIND FD_SPHINX_INSTALLS "man" indx)
  if (NOT indx EQUAL -1)
    install(
      DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/man
      OPTIONAL
      DESTINATION ${instdir}/man
      COMPONENT userdocs
      OPTIONAL
    )
  endif ()

endmacro()

