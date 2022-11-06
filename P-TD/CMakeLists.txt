cmake_minimum_required(VERSION 3.10)
project(P-TD)

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

add_definitions(-w) # omit warning information

add_executable(P-TD main.cpp index.cpp read.cpp update.cpp)

add_executable(P-TDP mainPost.cpp PostBoundary.cpp update.cpp index.cpp read.cpp)

if (Boost_FOUND)
    if(UNIX AND NOT APPLE)
        target_link_libraries(P-TD pthread boost_system boost_thread)# For linux
        target_link_libraries(P-TDP pthread boost_system boost_thread)# For linux
    endif()
    if(APPLE)
        target_link_libraries(P-TD boost_system boost_thread-mt)# For mac
        target_link_libraries(P-TDP boost_system boost_thread-mt)# For mac
    endif()
endif ()
