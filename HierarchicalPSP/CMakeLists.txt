cmake_minimum_required(VERSION 3.20)
project(Hierarchical_Index)

set(CMAKE_CXX_STANDARD 14)
# boost
set(BOOST_ROOT "/usr/lib/x86_64-linux-gnu")
find_package(Boost 1.71.0 COMPONENTS system filesystem REQUIRED)
if (Boost_FOUND)
    message(STATUS "Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
    message(STATUS "Boost_LIBRARIES: ${Boost_LIBRARY_DIRS}")
    message(STATUS "Boost_VERSION: ${Boost_VERSION}")
    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})
endif ()

if(UNIX AND NOT APPLE)
    include_directories(~/local/include)# For linux
    link_directories(~/local/lib)
    link_directories(/usr/lib/x86_64-linux-gnu)
endif()
if(APPLE)
    include_directories(/usr/include)# For mac
    link_directories(/usr/lib)
endif()
# add include path
##include_directories(/usr/include)
#include_directories(~/local/include)
## add library path
##link_directories(/usr/lib)
#link_directories(~/local/lib)
#link_directories(/usr/lib/x86_64-linux-gnu)

add_definitions(-w) # omit warning information

add_executable(Hier_DS main.cpp gtree.hpp gstartree.hpp)
#add_executable(Gtree gtree_build.cpp)

# link libraries
target_link_libraries(Hier_DS metis.a GKlib.a)
#target_link_libraries(Gtree metis.a GKlib.a)


if (Boost_FOUND)
    if(UNIX AND NOT APPLE)
        target_link_libraries(Hier_DS pthread boost_system boost_thread)# For linux
    endif()
    if(APPLE)
        target_link_libraries(Hier_DS boost_system boost_thread-mt)# For mac
    endif()
endif ()
