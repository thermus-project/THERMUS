# CMakeLists.txt in test

install(FILES all_predictions.C prediction.txt DESTINATION ${CMAKE_INSTALL_DOCDIR}/tests)
install(PROGRAMS inplace_test.sh DESTINATION ${CMAKE_INSTALL_DOCDIR}/tests)

add_subdirectory(expected_results)

