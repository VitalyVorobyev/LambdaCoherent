cmake_minimum_required (VERSION 3.1)
project (cohlam)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_BUILD_TYPE Release)

include_directories(
    include
    /home/vitaly/work/MyLibs/libLinAl/src
    /home/vitaly/Dist/podio/include
    /home/vitaly/Dist/fcc-edm/datamodel
)

link_directories(
    /home/vitaly/work/MyLibs/libLinAl/build
)

file(GLOB SOURCES src/*.cpp)

add_executable(${PROJECT_NAME} ${SOURCES})
target_link_libraries(${PROJECT_NAME}
    libpodio.so
    libdatamodel.so
    libLinal.so
)
