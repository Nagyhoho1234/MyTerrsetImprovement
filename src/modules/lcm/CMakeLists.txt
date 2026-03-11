add_library(aplaceholder_lcm STATIC
    ChangeAnalysis.cpp
    TransitionPotential.cpp
    MarkovChain.cpp
    CellularAutomata.cpp
    Prediction.cpp
    ChgAlloc.cpp
    LogisticReg.cpp
    Markov2.cpp
    FuzzyArtMap2.cpp
)

target_link_libraries(aplaceholder_lcm PUBLIC
    aplaceholder_core
    aplaceholder_io
)
