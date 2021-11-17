cmake_minimum_required(VERSION 3.14)
project(SSATAN-X_submit)

set(CMAKE_CXX_STANDARD 20)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")

FIND_PATH(LEMON_SOURCE_ROOT_DIR CMakeLists.txt
        PATHS ${CMAKE_SOURCE_DIR}/lib/lemon-1.3.1 ${CMAKE_SOURCE_DIR}/deps/lemon
        NO_DEFAULT_PATH
        DOC "Location of LEMON source as a CMAKE subproject")

IF(EXISTS ${LEMON_SOURCE_ROOT_DIR})
    ADD_SUBDIRECTORY(${LEMON_SOURCE_ROOT_DIR} deps/lemon)
    SET(LEMON_INCLUDE_DIRS
            ${LEMON_SOURCE_ROOT_DIR}
            ${CMAKE_BINARY_DIR}/deps/lemon
            )
    SET(LEMON_LIBRARIES lemon)
    UNSET(LEMON_ROOT_DIR CACHE)
    UNSET(LEMON_DIR CACHE)
    UNSET(LEMON_INCLUDE_DIR CACHE)
    UNSET(LEMON_LIBRARY CACHE)
ELSE()
    FIND_PACKAGE(LEMON QUIET NO_MODULE)
    FIND_PACKAGE(LEMON REQUIRED)
ENDIF()

## This line finds doxygen (for document creation)

FIND_PACKAGE(Doxygen)

INCLUDE_DIRECTORIES(
        ${PROJECT_SOURCE_DIR}
        ${PROJECT_BINARY_DIR}
        ${LEMON_INCLUDE_DIRS}
)

IF(CMAKE_COMPILER_IS_GNUCXX)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -static")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pthread")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

add_executable(SSATAN-X_submit main.cpp contact_network/Specie.cpp contact_network/Specie.h contact_network/ContactNetwork.cpp contact_network/ContactNetwork.h algorithms/SSA.cpp algorithms/SSA.h algorithms/SSATANX.cpp algorithms/SSATANX.h utilities/Utility.h utilities/Utility.cpp algorithms/AndersonTauLeap.h algorithms/AndersonTauLeap.cpp utilities/types.h nlohmann/json.h utilities/Settings.h utilities/Settings.cpp)