cmake_minimum_required (VERSION 3.8)

if (UNIX)
  add_definitions(-DLINUX)
endif (UNIX)
if (WIN32)
  add_definitions(-DWIN32)
endif (WIN32)

# Enable Hot Reload for MSVC compilers if supported.
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project ("blueseltzer")

# Add source to this project's executable.
add_executable (blueseltzer "main.cpp")

if (UNIX)
  find_package(JPEG REQUIRED)
  target_link_libraries(${PROJECT_NAME} JPEG::JPEG)
endif (UNIX)


if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET blueseltzer PROPERTY CXX_STANDARD 20)
endif()

# TODO: Add tests and install targets if needed.
