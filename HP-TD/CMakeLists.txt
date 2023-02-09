cmake_minimum_required(VERSION 3.10)
project(HP-TD)

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

if(UNIX AND NOT APPLE)
    include_directories(~/local/include)
    link_directories(~/local/lib)
endif()
if(APPLE)
    # add include path
    include_directories(/usr/include)
    # add library path
    link_directories(/usr/lib)
endif()


#link_directories(/usr/lib/x86_64-linux-gnu)

add_definitions(-w) # omit warning information

add_executable(HP-TD mainTD.cpp gtreeTD.hpp gstartreeTD.hpp H2H.hpp)
#add_executable(HP-DS main.cpp gtree.hpp gstartree.hpp)

# link libraries
target_link_libraries(HP-TD metis.a GKlib.a)
#target_link_libraries(HP-DS metis.a GKlib.a)

if (Boost_FOUND)
    if(UNIX AND NOT APPLE)
        target_link_libraries(HP-TD pthread boost_system boost_thread) #for linux
        #target_link_libraries(HP-DS pthread boost_system boost_thread) #for linux
    endif()
    if(APPLE)
        target_link_libraries(HP-TD boost_system boost_thread-mt) #for mac
        #target_link_libraries(HP-DS boost_system boost_thread-mt) #for mac
    endif()
endif ()
