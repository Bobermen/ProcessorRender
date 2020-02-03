cmake_minimum_required(VERSION 3.13)

project(ProcessorRendering LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(${CMAKE_SOURCE_DIR}/conancommand.cmake)
set_conan_params()

execute_process(COMMAND conan remote add bincrafters https://api.bintray.com/conan/bincrafters/public-conan)
if (WIN32)
    execute_process(COMMAND ${CMAKE_COMMAND} -E echo ${CONAN_INSTALL_COMMAND} COMMAND cmd)
elseif(UNIX)
    execute_process(COMMAND ${CMAKE_COMMAND} -E echo ${CONAN_INSTALL_COMMAND} COMMAND sh)
endif()

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)

conan_basic_setup(TARGET)


set(SOURCES
        main.cpp)

add_executable(ProcessorRendering ${SOURCES})

target_link_libraries(ProcessorRendering ${CONAN_LIBS} pthread)
