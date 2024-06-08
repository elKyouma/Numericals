file(GLOB SOURCES "*.cpp")

add_library(numericals ${SOURCES})
target_include_directories(numericals PUBLIC ${CMAKE_SOURCE_DIR}/include)
