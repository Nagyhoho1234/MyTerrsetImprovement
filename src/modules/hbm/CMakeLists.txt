add_library(aplaceholder_hbm STATIC
    SpeciesDistribution.cpp
    BiodiversityMetrics.cpp
    HabitatAssessment.cpp
    Habitat.cpp
    Jaccard.cpp
    OverlapUse.cpp
    SimWeight.cpp
)

target_link_libraries(aplaceholder_hbm PUBLIC
    aplaceholder_core
    aplaceholder_io
)
