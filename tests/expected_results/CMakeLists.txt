# CMakeLists.txt in test/expected_results

# test if ImageOS env variable is defined
if(DEFINED ENV{ImageOS})
    if($ENV{ImageOS} MATCHES macos14 OR $ENV{ImageOS} MATCHES macos15)
        install(FILES current_fit_lhc5020_fit_test_arm.txt  RENAME current_fit_lhc5020_fit_test.txt
            DESTINATION ${CMAKE_INSTALL_DOCDIR}/tests/expected_results)
    else()
        install(FILES current_fit_lhc5020_fit_test_x86_64.txt  RENAME current_fit_lhc5020_fit_test.txt
            DESTINATION ${CMAKE_INSTALL_DOCDIR}/tests/expected_results)
    endif()
else()
    install(FILES current_fit_lhc5020_fit_test_x86_64.txt  RENAME current_fit_lhc5020_fit_test.txt
        DESTINATION ${CMAKE_INSTALL_DOCDIR}/tests/expected_results)
endif()

install(FILES result.txt DESTINATION ${CMAKE_INSTALL_DOCDIR}/tests/expected_results)

