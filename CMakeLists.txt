CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)

PROJECT( TileGeneration )

SET( Boost_USE_STATIC_LIBS ON )
SET( Boost_NO_BOOST_CMAKE 1 )
FIND_PACKAGE( Boost 1.49 COMPONENTS system filesystem program_options REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIRS} )
LINK_DIRECTORIES( ${Boost_LIBRARY_DIRS} )

INCLUDE_DIRECTORIES("C:/Coding/Libraries/openslide-win32-20140125/include/openslide")
LINK_DIRECTORIES("C:/Coding/Libraries/openslide-win32-20140125/lib")

INCLUDE_DIRECTORIES("C:/Coding/Libraries/tiff-4.0.3/libtiff")
LINK_DIRECTORIES("C:/Coding/Libraries/tiff-4.0.3/libtiff")

SET( CurrentExe "TileGeneration" )
SET( CMAKE_C_FLAGS "-std=gnu99" )

SET( SOURCE_FILES
  TileGeneration.cpp
)

ADD_EXECUTABLE( ${CurrentExe} ${SOURCE_FILES})

TARGET_LINK_LIBRARIES( ${CurrentExe} ${Boost_LIBRARIES} libopenslide libtiff )
