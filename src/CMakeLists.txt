# add QuaCa library
set(quaca_sources
        Calculations/Integrations.cpp
        Friction/Friction.cpp
        GreensTensor/GreensTensor.cpp
        GreensTensor/GreensTensorFactory.cpp
        GreensTensor/GreensTensorPlate.cpp
        GreensTensor/GreensTensorPlateVacuum.cpp
        GreensTensor/GreensTensorVacuum.cpp
        Looper/Looper.cpp
        Looper/LooperFactory.cpp
        Looper/LooperV.cpp
        Looper/LooperZa.cpp
        MemoryKernel/MemoryKernelFactory.cpp
        MemoryKernel/OhmicMemoryKernel.cpp
		MemoryKernel/SinglePhononMemoryKernel.cpp
        Permittivity/PermittivityDrude.cpp
        Permittivity/PermittivityFactory.cpp
        Permittivity/PermittivityLorentz.cpp
        Polarizability/Polarizability.cpp
        PowerSpectrum/PowerSpectrum.cpp
        ReflectionCoefficients/ReflectionCoefficientsFactory.cpp
        ReflectionCoefficients/ReflectionCoefficientsLocBulk.cpp
        ReflectionCoefficients/ReflectionCoefficientsLocSlab.cpp
		version.cpp
        )
add_library(quaca SHARED ${quaca_sources})

# add git head and branch (taken from https://github.com/michael-hartmann/caps)
if (EXISTS "${PROJECT_SOURCE_DIR}/.git" AND ${GIT_FOUND})
	execute_process(
			COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
			WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
			OUTPUT_VARIABLE GIT_BRANCH
			OUTPUT_STRIP_TRAILING_WHITESPACE
	)

	execute_process(
			COMMAND ${GIT_EXECUTABLE} log -1 --format=%h
			WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
			OUTPUT_VARIABLE GIT_COMMIT_HASH
			OUTPUT_STRIP_TRAILING_WHITESPACE
	)

	target_compile_definitions(quaca PRIVATE -DGIT_HEAD="${GIT_COMMIT_HASH}" -DGIT_BRANCH="${GIT_BRANCH}")
	message("Git info added. commit hash: ${GIT_COMMIT_HASH} branch: ${GIT_BRANCH}")
endif ()

# link libraries
target_link_libraries(quaca
		${GSL_LIBRARY}
		${GSL_CBLAS_LIBRARY}
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
        ${Boost_FILESYSTEM_LIBRARY}
        ${Boost_SYSTEM_LIBRARY}
        ${Boost_PROGRAM_OPTIONS_LIBRARY}
        )

# set global header
set_target_properties(quaca PROPERTIES PUBLIC_HEADER include/Quaca.h)
