# ==============================================================================
# Project Configuration
# ==============================================================================
cmake_minimum_required (VERSION 3.8)

# Set project name
project ("MaxPR")

# Set C++17 standard (as specified in c_cpp_properties)
set (CMAKE_CXX_STANDARD 17)

# ==============================================================================
# Include Directories
# ==============================================================================
# Declare directory containing .h files
include_directories ("${CMAKE_SOURCE_DIR}/include")

# ==============================================================================
# Source Files
# ==============================================================================
# Find all .cpp files
file(GLOB SOURCE_FILES "${CMAKE_SOURCE_DIR}/src/*.cpp")

# ==============================================================================
# Build Target
# ==============================================================================
# Create executable
add_executable (main ${SOURCE_FILES})

# ==============================================================================
# Set working directory for Visual Studio (Debug & Release)
# ==============================================================================
set_target_properties(main PROPERTIES
    VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
)

# ==============================================================================
# Copy config and data files to build directory (auto-sync on every build)
# ==============================================================================
# Copy config folder
add_custom_command(TARGET main POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    "${CMAKE_SOURCE_DIR}/config"
    "${CMAKE_BINARY_DIR}/config"
    COMMENT "Syncing config files..."
)

# Copy data folder (only if not exists, to save time)
add_custom_command(TARGET main POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory_if_different
    "${CMAKE_SOURCE_DIR}/data"
    "${CMAKE_BINARY_DIR}/data"
    COMMENT "Syncing data files..."
)