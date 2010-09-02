########################################################

SET( LIBNAME kaltest)

########################################################


#------ ROOT stuff ----------
INCLUDE( ${PROJECT_SOURCE_DIR}/FindROOT.cmake )
SET( ENV{ROOTSYS} "${ROOT_HOME}" )
SET( ENV{PATH} $ENV{ROOTSYS}/bin:$ENV{PATH} )
SET( ENV{DYLD_LIBRARY_PATH} $ENV{ROOTSYS}/lib:$ENV{DYLD_LIBRARY_PATH} )
#-----------------------------

MESSAGE( STATUS " ROOT libraries:  ${ROOT_LIBRARIES} " )




################ external ilcsoft packages ############################
INCLUDE( "MacroLoadPackage" )

#LOAD_PACKAGE( KalTest REQUIRED )
LOAD_PACKAGE( GEAR REQUIRED )
LOAD_PACKAGE( Marlin REQUIRED )
IF( APPLE )
LOAD_PACKAGE( streamlog REQUIRED )
ENDIF( APPLE ) 


###### include directories #######################################

#### FIXME: need KalTest  finder module ########################
INCLUDE_DIRECTORIES( ${ROOT_INCLUDE_DIRS} )

IF( NOT DEFINED KalTest_HOME )
MESSAGE( FATAL_ERROR "KalTest_HOME not set : please set it with cmake -D KalTest_HOME=...")
ENDIF()
INCLUDE_DIRECTORIES( "${KalTest_HOME}/include" )

INCLUDE_DIRECTORIES( . )

#INCLUDE_DIRECTORIES( ${Marlin_INCLUDE_DIR} )
#INCLUDE_DIRECTORIES( ${GEAR_INCLUDE_DIR} )

#####################################################################

aux_source_directory( . lib_srcs)

ADD_LIBRARY( ${LIBNAME}  ${lib_srcs} )

TARGET_LINK_LIBRARIES( ${LIBNAME} ${KalTest_HOME}/lib/libS4Kalman.dylib)
TARGET_LINK_LIBRARIES( ${LIBNAME} ${KalTest_HOME}/lib/libS4Kaltrack.dylib)
TARGET_LINK_LIBRARIES( ${LIBNAME} ${KalTest_HOME}/lib/libS4KalGeom.dylib)
TARGET_LINK_LIBRARIES( ${LIBNAME} ${KalTest_HOME}/lib/libS4Utils.dylib)

TARGET_LINK_LIBRARIES( ${LIBNAME} ${ROOT_LIBRARIES})
#TARGET_LINK_LIBRARIES( ${LIBNAME} ${GEAR_LIBRARIES})
#TARGET_LINK_LIBRARIES( ${LIBNAME} ${Marlin_LIBRARIES})
#TARGET_LINK_LIBRARIES( ${LIBNAME} ${streamlog_LIBRARIES})

#####################################################################

INSTALL( TARGETS ${LIBNAME} DESTINATION lib PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE )

MESSAGE( STATUS " --- kaltest headers :  ${PROJECT_SOURCE_DIR}/kaltest/EXKalTest.h" )

INSTALL( FILES "${PROJECT_SOURCE_DIR}/kaltest/KalTest.h" DESTINATION ${PROJECT_SOURCE_DIR}/include )
