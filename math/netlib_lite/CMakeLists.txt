

# SEK:  Will use lsode later for SOL current idea.  
add_subdirectory(lsode)
#add_subdirectory(other)
add_subdirectory(r8slatec)

# This is not working for some reason
set(math_libs lsode r8slatec ${math_libs})
