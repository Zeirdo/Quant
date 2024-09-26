# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Debug")
  file(REMOVE_RECURSE
  "CMakeFiles/topmg_autogen.dir/AutogenUsed.txt"
  "CMakeFiles/topmg_autogen.dir/ParseCache.txt"
  "topmg_autogen"
  )
endif()
