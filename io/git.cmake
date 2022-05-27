######################################################################
#
# Parse information about the git repo.
#
######################################################################

# Get the current working branch
execute_process(
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the latest abbreviated commit hash of the working branch
execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the origin url
execute_process(
  COMMAND git config --get remote.origin.url
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_URL
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Check for modifications
execute_process(
  COMMAND bash -c "git diff-index --quiet HEAD || echo 'with local modifications'"
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_MODIFICATIONS
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get date of the working branch
execute_process(
  COMMAND git log -1 --format=%ai
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_DATE
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the date in a different format
execute_process(
  COMMAND git log -1 --format=%cD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_DATE_ALT
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
