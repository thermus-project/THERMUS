# CMakeLists.txt in test

install(FILES prediction.C prediction.txt prediction_expected.txt DESTINATION ${CMAKE_INSTALL_DOCDIR}/tests)
install(PROGRAMS inplace_test.sh DESTINATION ${CMAKE_INSTALL_DOCDIR}/tests)

