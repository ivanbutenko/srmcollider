project("SRMCollider")
cmake_minimum_required(VERSION 2.6)

set(CMAKE_BUILD_TYPE Release)

# All executables 
set(my_executables
testPrecursor
testgetNonUIS
testIntegrated
testRangetree
testeUIS
testLibrary
combinatorics_test
run_integrated
)

# All classes that will go into the library
set(my_sources
  calculate_eUIS.cpp
  combinatorics.h
  srmcolliderLib.cpp
  integratedrun.cpp
  rangetree.cpp
)

# CGAL is used header-only, no need for the libraries
find_package(CGAL REQUIRED)
message("Include dirs of CGAL: " ${CGAL_INCLUDE_DIRS} )
# message("Libs of CGAL: " ${CGAL_LIBRARIES} )

find_package(PythonLibs REQUIRED)
message("Include dirs of Python: " ${PYTHON_INCLUDE_DIRS} )
message("Libs of Python: " ${PYTHON_LIBRARIES} )

find_package(Boost COMPONENTS filesystem system unit_test_framework python REQUIRED)
message("Include dirs of boost: " ${Boost_INCLUDE_DIRS} )
message("Libs of boost: " ${Boost_LIBRARIES} )

include_directories(
  ${Boost_INCLUDE_DIRS}
  ${PYTHON_INCLUDE_DIRS}
  ${CGAL_INCLUDE_DIRS}
)

# cgal requires the -frounding-math option
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -frounding-math -Wall") 

add_library(SRMCollider SHARED ${my_sources})

## add targets for the executables
foreach(i ${my_executables})
  add_executable(${i} ${i}.C)
  target_link_libraries(${i} SRMCollider ${PYTHON_LIBRARIES} ${Boost_LIBRARIES} )
endforeach(i)

ENABLE_TESTING()
ADD_TEST(testgetNonUIS testgetNonUIS)
ADD_TEST(testPrecursor testPrecursor)
ADD_TEST(testIntegrated testIntegrated)
ADD_TEST(testRangetree testRangetree)
ADD_TEST(testeUIS testeUIS)
ADD_TEST(testLibrary testLibrary)
ADD_TEST(combinatorics_test combinatorics_test)

