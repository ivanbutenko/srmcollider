project("SRMCollider")
cmake_minimum_required(VERSION 2.6)

# set(CMAKE_BUILD_TYPE Release)

## list all your executables here (a corresponding .C file should exist, e.g. Main.C)
set(my_executables
testmain
testIntegrated
)


#find_package(CGAL)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}") 
include_directories("/usr/include/python2.6")

## add targets for the executables
foreach(i ${my_executables})
  add_executable(${i} ${i}.C)
  ## link executables against OpenMS and also OpenSWATH
  target_link_libraries(${i} boost_python boost_unit_test_framework python2.6)
endforeach(i)


ENABLE_TESTING()
ADD_TEST(testmain testmain)
ADD_TEST(testIntegrated testIntegrated)

