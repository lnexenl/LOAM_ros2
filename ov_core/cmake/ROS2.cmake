cmake_minimum_required(VERSION 3.3)

# Find ros dependencies
find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(cv_bridge REQUIRED)

# Describe ROS project
option(ENABLE_ROS "Enable or disable building with ROS (if it is found)" ON)
if (NOT ENABLE_ROS)
    message(FATAL_ERROR "Build with ROS1.cmake if you don't have ROS.")
endif ()
add_definitions(-DROS_AVAILABLE=2)

# Include our header files
include_directories(
        src
        ${EIGEN3_INCLUDE_DIR}
        ${Boost_INCLUDE_DIRS}
)

# Set link libraries used by all binaries
list(APPEND thirdparty_libraries
        ${Boost_LIBRARIES}
        ${OpenCV_LIBRARIES}
        ${spdlog_LIBRARIES}
)

##################################################
# Make the core library
##################################################

list(APPEND LIBRARY_SOURCES
        src/dummy.cpp
        src/cpi/CpiV1.cpp
        src/cpi/CpiV2.cpp
        src/types/Landmark.cpp
        src/utils/print.cpp
        src/state/State.cpp
        src/state/StateHelper.cpp
)
file(GLOB_RECURSE LIBRARY_HEADERS "src/*.h")
add_library(ov_core_lib SHARED ${LIBRARY_SOURCES} ${LIBRARY_HEADERS})
ament_target_dependencies(ov_core_lib rclcpp cv_bridge)
target_link_libraries(ov_core_lib ${thirdparty_libraries})
target_include_directories(ov_core_lib PUBLIC src/)
install(TARGETS ov_core_lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin
        PUBLIC_HEADER DESTINATION include
)
install(DIRECTORY src/
        DESTINATION include
        FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp"
)
ament_export_include_directories(include)
ament_export_libraries(ov_core_lib)

# finally define this as the package
ament_package()