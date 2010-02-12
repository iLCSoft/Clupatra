###############################################
# cmake configuration file for clupatra
# @author Jan Engels, DESY
###############################################

SET( clupatra_FOUND FALSE )
MARK_AS_ADVANCED( clupatra_FOUND )

# do not store find results in cache
SET( clupatra_INCLUDE_DIR clupatra_INCLUDE_DIR-NOTFOUND )

FIND_PATH( clupatra_INCLUDE_DIR
	NAMES clupatra.h
	PATHS /Users/gaede/marlin/Clupatra
	PATH_SUFFIXES include
	NO_DEFAULT_PATH
)
IF( NOT clupatra_INCLUDE_DIR )
    MESSAGE( STATUS "Check for clupatra: ${clupatra_HOME}"
					" -- failed to find clupatra include directory!!" )
ELSE( NOT clupatra_INCLUDE_DIR )
    MARK_AS_ADVANCED( clupatra_INCLUDE_DIR )
ENDIF( NOT clupatra_INCLUDE_DIR )


# do not store find results in cache
SET( clupatra_LIB clupatra_LIB-NOTFOUND )

FIND_LIBRARY( clupatra_LIB
	NAMES clupatra
	PATHS /Users/gaede/marlin/Clupatra
	PATH_SUFFIXES lib
	NO_DEFAULT_PATH
)
IF( NOT clupatra_LIB )
    MESSAGE( STATUS "Check for clupatra: ${clupatra_HOME}"
					" -- failed to find clupatra library!!" )
ELSE( NOT clupatra_LIB )
    MARK_AS_ADVANCED( clupatra_LIB )
ENDIF( NOT clupatra_LIB )


# set variables and display results
IF( clupatra_INCLUDE_DIR AND clupatra_LIB )
    SET( clupatra_FOUND TRUE )
    SET( clupatra_INCLUDE_DIRS ${clupatra_INCLUDE_DIR} )
    SET( clupatra_LIBRARY_DIRS "/Users/gaede/marlin/Clupatra/lib" )
	SET( clupatra_LIBRARIES ${clupatra_LIB} )
    MARK_AS_ADVANCED( clupatra_INCLUDE_DIRS clupatra_LIBRARY_DIRS clupatra_LIBRARIES )
	MESSAGE( STATUS "Check for clupatra: ${clupatra_HOME} -- works" )
ELSE( clupatra_INCLUDE_DIR AND clupatra_LIB )
	IF( clupatra_FIND_REQUIRED )
		MESSAGE( FATAL_ERROR "Check for clupatra: ${clupatra_HOME} -- failed!!" )
    ELSE( clupatra_FIND_REQUIRED )
        MESSAGE( STATUS "Check for clupatra: ${clupatra_HOME}"
						" -- failed!! will skip this package..." )
    ENDIF( clupatra_FIND_REQUIRED )
ENDIF( clupatra_INCLUDE_DIR AND clupatra_LIB )
