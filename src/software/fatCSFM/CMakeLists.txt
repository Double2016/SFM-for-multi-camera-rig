
include_directories(CMAKE_SOURCE_DIR)

find_package(SIFTGPU REQUIRED)
set(SIFTGPU_LIBRARY ${SIFTGPU_DIR}/lib/siftgpu.lib)
set(SIFTGPU_INCLUDE ${SIFTGPU_DIR}/include)
#find_package(BOOST REQUIRED)
find_package(CGAL REQUIRED)
include(${CGAL_USE_FILE})
#set(CGAL_INCLUDE ${CGAL_DIR}/include)
#set(CGAL_LIBRARY ${CGAL_DIR}/lib/CGAL-vc140-mt-4.13.lib)


find_package(OpenGL REQUIRED)
find_package(Glew REQUIRED)

message("asdasdwqeqweqwezd************************************** ${CGAL_LIBRARIES}")
#${SIFTGPU_LIBRARY} ${OPENGL_LIBRARIES} ${GLEW_LIBRARIES}
ADD_EXECUTABLE(openMVG_main_fatCSFM    main_fatCSFM.cpp fatCSFM.hpp siftgpu_GLSL.h)
target_link_libraries(openMVG_main_fatCSFM
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
	${OPENGL_LIBRARIES}
	${GLEW_LIBRARIES}
	${SIFTGPU_LIBRARY}
	${CGAL_LIBRARIES}
)

target_include_directories(openMVG_main_fatCSFM
PRIVATE
	${SIFTGPU_INCLUDE}
	${CGAL_INCLUDE}
	${Boost_INCLUDE_DIR}
)
# Installation rules
if (MSVC)
  set_target_properties(openMVG_main_fatCSFM PROPERTIES COMPILE_FLAGS "/bigobj")
endif (MSVC)
set_property(TARGET openMVG_main_fatCSFM PROPERTY FOLDER OpenMVG/software)
install(TARGETS openMVG_main_fatCSFM DESTINATION bin/)

