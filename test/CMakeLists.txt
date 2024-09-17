# CMakeLists.txt in test

install(FILES prediction.txt prediction_expected.txt DESTINATION ${CMAKE_INSTALL_DOCDIR}/tests)
install(FILES prediction_dest.C DESTINATION ${CMAKE_INSTALL_DOCDIR}/tests RENAME prediction.C)

