add_executable(matrix_solver_tests matrix_solver_tests.cpp)
target_include_directories(matrix_solver_tests PUBLIC ${CMAKE_SOURCE_DIR}/include)

target_link_libraries(matrix_solver_tests
 PRIVATE
  GTest::GTest
  matrix_solver)

add_test(matrix_gtests matrix_solver_tests_test)
