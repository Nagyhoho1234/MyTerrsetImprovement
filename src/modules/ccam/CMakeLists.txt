add_library(aplaceholder_ccam STATIC
    ClimateAdaptation.cpp
    EcoCrop.cpp
    ScenGen.cpp
)

target_link_libraries(aplaceholder_ccam PUBLIC
    aplaceholder_core
    aplaceholder_io
)
