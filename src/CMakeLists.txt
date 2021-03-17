cmake_minimum_required (VERSION 3.8)

set (mason2_args ${APP_TEMPLATE_EXTERNAL_PROJECT_CMAKE_ARGS})
list (APPEND mason2_args "-DSEQAN_BUILD_SYSTEM=APP:mason2")

include (ExternalProject)
ExternalProject_Add (
    mason2
    PREFIX mason2
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/lib/seqan"
    CMAKE_ARGS "${mason2_args}"
    INSTALL_COMMAND ""
)
unset (mason2_args)

ExternalProject_Get_property(mason2 BINARY_DIR)
install(DIRECTORY "${BINARY_DIR}/bin/" DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}" USE_SOURCE_PERMISSIONS)

add_library (common_interface INTERFACE)
target_link_libraries ("common_interface" INTERFACE seqan3::seqan3)
target_include_directories ("common_interface" INTERFACE ../include)

add_executable ("generate_reads" generate_reads.cpp)
target_link_libraries ("generate_reads" "common_interface")

add_executable ("split_sequence" split_sequence.cpp)
target_link_libraries ("split_sequence" "common_interface")

add_executable ("count_minimiser" count_minimiser.cpp)
target_link_libraries ("count_minimiser" "common_interface")

add_executable ("apply_taxsbp" apply_taxsbp.cpp)
target_link_libraries ("apply_taxsbp" "common_interface")

add_executable ("generate_reads_refseq" generate_reads_refseq.cpp)
target_link_libraries ("generate_reads_refseq" "common_interface")
