add_library(aplaceholder_esm STATIC
    EcosystemServices.cpp
    Carbon.cpp
    Pollination.cpp
    CoastalVuln.cpp
    WaveEnergy.cpp
    Sediment.cpp
    SoilSalinity.cpp
)

target_link_libraries(aplaceholder_esm PUBLIC
    aplaceholder_core
    aplaceholder_io
)
