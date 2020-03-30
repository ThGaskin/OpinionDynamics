# Add the model target
add_model(OpDyn OpDyn.cc)
# NOTE The target should have the same name as the model folder and the *.cc
# Add test directories
add_subdirectory(tests EXCLUDE_FROM_ALL)

