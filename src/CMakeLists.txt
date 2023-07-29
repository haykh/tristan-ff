# if userfile not defined
# list all files in user directory
set(user "" CACHE STRING "userfile")
if(user STREQUAL "")
  file(GLOB_RECURSE userfiles ${CMAKE_CURRENT_SOURCE_DIR}/../user/*.F90)
  message(STATUS "Please specify userfile:")
  foreach(userfile ${userfiles})
    get_filename_component(userfile ${userfile} NAME)
    string(REGEX REPLACE ".F90$" "" userfile ${userfile})
    message(STATUS ": ${userfile}")
  endforeach()
  message(FATAL_ERROR "userfile not defined")
else()
  set(userfile "${CMAKE_CURRENT_SOURCE_DIR}/../user/${user}.F90")

  set(executable ${PROJECT_NAME}.xc)
  file(GLOB_RECURSE SOURCES *.F90)
  add_executable(${executable} ${SOURCES} ${userfile})
  add_definitions(-D NGHOST=2)

  target_link_libraries(${executable} HDF5::HDF5)
endif()
