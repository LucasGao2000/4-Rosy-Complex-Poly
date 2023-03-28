cmake_minimum_required(VERSION 3.23)
project(4-ROSY VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 20)

include(CTest)
enable_testing()

SET(OPENMESH_PATH "/usr/local/include")
SET(OPENMESH_LIB_PATH "/usr/local/lib")
set(EIGEN "extern/Eigen")

include_directories(${PROJECT_SOURCE_DIR}/extern)
include_directories(${OPENMESH_PATH})

link_directories(${OPENMESH_LIB_PATH})

add_executable(4-ROSY 
    main.cpp
    src/FrameField.cpp
)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

target_link_libraries(${PROJECT_NAME}
        OpenMeshCore
        OpenMeshTools
        )

target_compile_definitions(${PROJECT_NAME} PUBLIC -DIS64BITPLATFORM)

if(${CMAKE_BUILD_TYPE} MATCHES "Release")
    message(STATUS "Release版本")
    set(BuildType "Release")
    set(CMAKE_C_FLAGS "$ENV{CFLAGS} -DNODEBUG -O3 -Wall")
    set(CMAKE_CXX_FLAGS "$ENV{CXXFLAGS} -DNODEBUG -O3 -Wall")
else()
    message(STATUS "Debug版本")
    set(BuildType "Debug")
    set(CMAKE_C_FLAGS "$ENV{CFLAGS} -O0 -g -Wall")
    set(CMAKE_CXX_FLAGS "$ENV{CXXFLAGS} -O0 -g -Wall")
endif()
