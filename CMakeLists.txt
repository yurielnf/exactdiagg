cmake_minimum_required(VERSION 3.9)

project(exactdiagg LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


set(SOURCES 
    main.cpp
    examples/example_hernan/example_hernan.cpp
    examples/example_nair/example_cadenaAA3.cpp
    examples/example_nair/example_cadenaAA5Open.cpp
    examples/parameters.cpp
    examples/example_tb.cpp
    examples/example_hubbard.cpp
    examples/example_hubbard2.cpp
    examples/example_2qd.cpp
)

set(HEADERS
    exactdiagg/hamsolver.h
    exactdiagg/lanczos.h
    exactdiagg/qoperator.h
    exactdiagg/fockbasis.h
    exactdiagg/symmetrygroup.h
    exactdiagg/all.h
    examples/all.h
    examples/example_nair/cadenitaaa3.h
    examples/example_nair/cadenitaaa5open.h
    examples/parameters.h
)

add_executable(${PROJECT_NAME} 
    ${HEADERS}
    ${SOURCES}
)

target_include_directories(${PROJECT_NAME} PRIVATE .)

#set_property(TARGET ${PROJECT_NAME} PROPERTY POSITION_INDEPENDENT_CODE ON)


#---------- external dependencies ---------

find_package(Armadillo REQUIRED)
target_include_directories(${PROJECT_NAME} PRIVATE ${ARMADILLO_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PRIVATE ${ARMADILLO_LIBRARIES})
