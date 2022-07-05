# Author: Hyunbin Kim (khb7840@gmail.com)
set (CMAKE_CXX_STANDARD 17)
# 01. Collect files
file(GLOB_RECURSE SOURCE_CPP ${CMAKE_CURRENT_SOURCE_DIR}/ *.cpp)
file(GLOB_RECURSE SOURCE_H ${CMAKE_CURRENT_SOURCE_DIR}/ *.h)
# 02. Add executable
add_library(${PROJECT_NAME} OBJECT ${SOURCE_CPP} ${SOURCE_H})