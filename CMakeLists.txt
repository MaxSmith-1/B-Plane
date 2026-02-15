# Configure Project
cmake_minimum_required(VERSION 3.10)
project(B_Plane VERSION 1.0.0 LANGUAGES CXX)

set(BUILD_TESTING OFF CACHE BOOL "Disable testing" FORCE)

# Set default build type if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find packages - Use Conan's target names
find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED)
find_package(jsoncpp REQUIRED)
find_package(CURL REQUIRED)


# Generate executable
add_executable(B_Plane 
    main.cpp
    src/Simulator.cpp
)

# Link libraries - Use lowercase target names for Conan packages
target_link_libraries(B_Plane 
    PRIVATE 
    Eigen3::Eigen
    boost::boost
    JsonCpp::JsonCpp
    CURL::libcurl
)

target_include_directories(B_Plane PRIVATE src)

# Copy executable
if(WIN32)
    add_custom_command(TARGET B_Plane POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy 
            $<TARGET_FILE:B_Plane>
            ${CMAKE_SOURCE_DIR}/executable/B_Plane.exe
            ${CMAKE_SOURCE_DIR}/executable/
        COMMAND ${CMAKE_COMMAND} -E remove 
            $<TARGET_FILE:B_Plane>
    )
else()
    # For Linux
    add_custom_command(TARGET B_Plane POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy 
            $<TARGET_FILE:B_Plane>
            ${CMAKE_SOURCE_DIR}/executable/B_Plane
            ${CMAKE_SOURCE_DIR}/executable/
        COMMAND ${CMAKE_COMMAND} -E remove 
            $<TARGET_FILE:B_Plane>
    )
endif()