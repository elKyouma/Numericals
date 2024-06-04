add_executable(numericals_tests numericals_tests.cpp)
target_include_directories(numericals_tests PUBLIC ${CMAKE_SOURCE_DIR}/include)

target_link_libraries(numericals_tests
 PRIVATE
  GTest::GTest
  numericals)

add_test(numericals_gtests numericals_test)
