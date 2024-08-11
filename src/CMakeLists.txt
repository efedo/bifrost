file(GLOB sources *.cpp roaring.c)
file(GLOB headers *.h *.hpp *.hh *.tcc)

list(REMOVE_ITEM sources Bifrost.cpp)

add_definitions(-DMAX_KMER_SIZE=${MAX_KMER_SIZE})
add_definitions(-DMAX_GMER_SIZE=${MAX_GMER_SIZE})

add_library(bifrost_static STATIC ${sources} ${headers})
add_library(bifrost_dynamic SHARED ${sources} ${headers})

set_target_properties(bifrost_static PROPERTIES OUTPUT_NAME "bifrost_static")
set_target_properties(bifrost_dynamic PROPERTIES OUTPUT_NAME "bifrost_dynamic")

target_include_directories(bifrost_static PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(bifrost_dynamic PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(Bifrost Bifrost.cpp)

find_package(Threads REQUIRED)
target_link_libraries(bifrost_static ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(bifrost_dynamic ${CMAKE_THREAD_LIBS_INIT})

#find_package(ZLIB REQUIRED)

include(FetchContent)

set(ZLIB_COMPAT ON)
set(ZLIB_ENABLE_TESTS OFF)

FetchContent_Declare(
  zlib-ng
  GIT_REPOSITORY https://github.com/zlib-ng/zlib-ng.git
  GIT_TAG        d54e3769be0c522015b784eca2af258b1c026107 # 2.2.1
  OVERRIDE_FIND_PACKAGE
)
FetchContent_MakeAvailable(zlib-ng)


target_link_libraries(bifrost_static zlib-ng)
target_link_libraries(bifrost_dynamic zlib-ng)

#if (ZLIB_FOUND)
#    include_directories(${ZLIB_INCLUDE_DIRS})
#else()
#    message(FATAL_ERROR "zlib not found. Required for to output files")
#endif(ZLIB_FOUND)

target_link_libraries(Bifrost bifrost_dynamic)

install(TARGETS Bifrost DESTINATION bin)
install(TARGETS bifrost_dynamic DESTINATION lib)
install(TARGETS bifrost_static DESTINATION lib)
install(FILES ${headers} DESTINATION include/bifrost)
