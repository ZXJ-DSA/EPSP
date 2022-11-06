cmake_minimum_required(VERSION 3.10)
project(P-PLL)

set(CMAKE_CXX_STANDARD 11)


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

add_executable(P-PLL main.cpp Function.cpp Update.cpp PSL.cpp)
add_executable(P-PLLP mainPost.cpp Function.cpp Update.cpp PSL.cpp PostBoundary.cpp)

if (Boost_FOUND)
    if(UNIX AND NOT APPLE)
        target_link_libraries(P-PLL pthread boost_system boost_thread)# For linux
        target_link_libraries(P-PLLP pthread boost_system boost_thread)# For linux
    endif()
    if(APPLE)
        target_link_libraries(P-PLL boost_system boost_thread-mt)# For mac
        target_link_libraries(P-PLLP boost_system boost_thread-mt)# For mac
    endif()
endif ()
