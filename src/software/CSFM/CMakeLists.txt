
include_directories(CMAKE_SOURCE_DIR)


ADD_EXECUTABLE(openMVG_main_CSFM    main_CSFM.cpp CSFM.hpp)
target_link_libraries(openMVG_main_CSFM
  PRIVATE
    openMVG_system
    openMVG_image
    openMVG_features
    openMVG_sfm
    openMVG_exif
    ${STLPLUS_LIBRARY}
	openMVG_multiview
	vlsift
	openMVG_graph
	openMVG_matching_image_collection
)
#target_include_directories(openMVG_main_CSFM
#)
# Installation rules
if (MSVC)
  set_target_properties(openMVG_main_CSFM PROPERTIES COMPILE_FLAGS "/bigobj")
endif (MSVC)
set_property(TARGET openMVG_main_CSFM PROPERTY FOLDER OpenMVG/software)
install(TARGETS openMVG_main_CSFM DESTINATION bin/)

