add_model_tests(# Use consistent capitalization for the model name!
                MODEL_NAME OpDyn
                # The sources of the model tests to carry out. Each of these
                # will become a test target with the same name.
                SOURCES
                    "test_utils.cc"
                    "test_ageing.cc"
                # Optional: Files to be copied to the build directory
                AUX_FILES
                    "test_config.yml"
                )
