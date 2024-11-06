# CMakeLists.txt in test

install(FILES all_predictions.C prediction.txt lhc5020_fit_charm.C lhc5020_final_0_single_charm.txt DESTINATION ${CMAKE_INSTALL_DOCDIR}/tests)
install(PROGRAMS inplace_test.sh DESTINATION ${CMAKE_INSTALL_DOCDIR}/tests)

add_subdirectory(expected_results)

