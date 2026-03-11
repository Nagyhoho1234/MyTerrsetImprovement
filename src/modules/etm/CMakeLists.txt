add_library(aplaceholder_etm STATIC
    TrendAnalysis.cpp
    SeasonalDecomp.cpp
    TimeSeriesStats.cpp
    Cca.cpp
    Eot.cpp
    Hants.cpp
    Modwt.cpp
    Mssa.cpp
    Sta.cpp
    MultiReg.cpp
    Mlr.cpp
    BioClimVars.cpp
    Deseason.cpp
    TAggregate.cpp
    TDiff.cpp
    TCor.cpp
    HarmonicReg.cpp
)

target_link_libraries(aplaceholder_etm PUBLIC
    aplaceholder_core
    aplaceholder_io
)
