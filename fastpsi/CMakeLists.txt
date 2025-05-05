

add_executable(fastpsi  fastpsi.cpp util.cpp)

target_link_libraries(fastpsi volePSI_Tests)

if(VOLE_PSI_NO_LTO)
    target_link_options(fastpsi PUBLIC "-fno-lto")
endif()