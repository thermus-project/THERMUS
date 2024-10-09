# Bash and Root scripts

# Convert run_thermus.in to run_thermus
set(THERMUS ${CMAKE_INSTALL_PREFIX})
set(THERMUS_LIB ${CMAKE_INSTALL_LIBDIR})
configure_file(run_thermus.in run_thermus)

# Install run_thermus
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/run_thermus DESTINATION ${CMAKE_INSTALL_BINDIR})

# Install use_thermus.C
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/use_thermus.C DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${CMAKE_PROJECT_NAME}/root/)
