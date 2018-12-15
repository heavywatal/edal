file(GLOB srcfiles *.cpp)
add_library(objlib OBJECT ${srcfiles})
target_include_directories(objlib INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
)
target_link_libraries(objlib PUBLIC
  wtl::wtl
  wtl::sfmt
  clippson::clippson
  ZLIB::ZLIB
)