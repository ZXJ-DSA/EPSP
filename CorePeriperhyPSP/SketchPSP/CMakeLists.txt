cmake_minimum_required(VERSION 3.10)
project(BatchHL)

set(CMAKE_CXX_STANDARD 14)

# boost
find_package(Boost 1.71.0 COMPONENTS system filesystem REQUIRED)
if (Boost_FOUND)
    message(STATUS "Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
    message(STATUS "Boost_LIBRARIES: ${Boost_LIBRARY_DIRS}")
    message(STATUS "Boost_VERSION: ${Boost_VERSION}")
    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})
endif ()

add_definitions(-w)

#add_executable(BatchHL Code/main.cpp Code/BatchHL.h Code/two_layer_queue.h)
#add_executable(BatchHLW Code/mainW.cpp Code/BatchHLW.hpp Code/BatchHLW.h)
add_executable(BHL mainWW.cpp BatchHLWW.hpp BatchHLWW.h)
if (Boost_FOUND)
    target_link_libraries(BHL pthread boost_system boost_thread)# For linux
	#target_link_libraries(BHL boost_system boost_thread-mt)# For mac
endif ()

