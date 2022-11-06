cmake_minimum_required(VERSION 3.10)
project(HP-DS)

set(CMAKE_CXX_STANDARD 14)
# boost
#set(BOOST_ROOT "/usr/lib/x86_64-linux-gnu")
find_package(Boost 1.71.0 COMPONENTS system filesystem REQUIRED)
if (Boost_FOUND)
    message(STATUS "Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
    message(STATUS "Boost_LIBRARIES: ${Boost_LIBRARY_DIRS}")
    message(STATUS "Boost_VERSION: ${Boost_VERSION}")
    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})
endif ()

# add include path
include_directories(/usr/include)
#include_directories(~/local/include)

# add library path
link_directories(/usr/lib)
#link_directories(~/local/lib)
#link_directories(/usr/lib/x86_64-linux-gnu)

add_definitions(-w) # omit warning information

add_executable(HP-DS main.cpp gtree.hpp gstartree.hpp)

# link libraries
target_link_libraries(HP-DS metis.a GKlib.a)

if (Boost_FOUND)
    if(UNIX AND NOT APPLE)
        target_link_libraries(HP-DS pthread boost_system boost_thread) #for linux
    endif()
    if(APPLE)
        target_link_libraries(HP-DS boost_system boost_thread-mt) #for mac
    endif()
endif ()
