add_library(NPtools SHARED)

target_sources(NPtools PRIVATE
    "RandomNumber.cpp"
    "InputParser.cpp"
    "CommandLineArguments.cpp"
    "FileSystem.cpp"
)

set(NPTOOLS_HEADER
    "Assert.h"
    "CommonTypes.h"
    "GenericFactory.h"
    "InputParser.h"
    "RandomNumber.h"
    "CommandLineArguments.h"
    "Registry.h"
    "FileSystem.h"
)

target_include_directories(NPtools
                        PUBLIC 
                        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
                        $<INSTALL_INTERFACE:${PROJECT_BINARY_DIR}/include>
)
# target_compile_options(NPtools PRIVATE -fPIC)


include(InstallLibraryWithStandardSetup)
InstallLibraryWithStandardSetup(NPtools 
    INCLUDE_SUBDIR NPtools
    INCLUDES ${NPTOOLS_HEADER})

