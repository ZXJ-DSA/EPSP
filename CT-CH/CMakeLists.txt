cmake_minimum_required(VERSION 3.10)
project(CT3)

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

#add_executable(CT-DS mainCTDijkstra.cpp CTpartition.cpp BasicFun.cpp PSL3.cpp PeripheryDijk.cpp)
#add_executable(CT-TD mainCTH2H.cpp CTpartition.cpp BasicFun.cpp PSL.cpp PeripheryH2H.cpp H2H.cpp)
#add_executable(CT-TD2 mainCTH2H2.cpp CTpartition2.cpp BasicFun2.cpp PSL2.cpp PeripheryH2H2.cpp H2H2.cpp)#simplified version
#add_executable(CT-TD3 mainCTH2H3.cpp CTpartition3.cpp BasicFun3.cpp PSL3.cpp PLL3.cpp PeripheryH2H3.cpp H2H3.cpp)#simplified version
add_executable(CT-CH mainCTCH.cpp CTpartition3.cpp BasicFun3.cpp PSL3.cpp PeripheryH2H3.cpp H2H3.cpp CH.cpp)#simplified version

if (Boost_FOUND)
    if(UNIX AND NOT APPLE)
#        target_link_libraries(CT-DS pthread boost_system boost_thread)# For linux
#        target_link_libraries(CT-TD pthread boost_system boost_thread)# For linux
#        target_link_libraries(CT-TD2 pthread boost_system boost_thread)# For linux
#        target_link_libraries(CT-TD3 pthread boost_system boost_thread)# For linux
	    target_link_libraries(CT-CH pthread boost_system boost_thread)# For linux
    endif()
    if(APPLE)
#        target_link_libraries(CT-DS boost_system boost_thread-mt)# For mac
#        target_link_libraries(CT-TD boost_system boost_thread-mt)# For mac
#        target_link_libraries(CT-TD2 boost_system boost_thread-mt)# For mac
#        target_link_libraries(CT-TD3 boost_system boost_thread-mt)# For mac
	    target_link_libraries(CT-CH boost_system boost_thread-mt)# For mac
    endif()
endif ()
