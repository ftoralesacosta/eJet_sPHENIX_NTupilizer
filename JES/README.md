JES_myjetanalysis.cc runs the JES/JER analysis on TTrees created from myjetanalysis from the sPHENIX tutorials repository
JES_TreeMaker run the JES/JER analysis on TTrees created from https://github.com/sPHENIX-Collaboration/macros/pull/222

The converter code simply take trees from TreeMaker output, and converts them to TTree similar to myjetanalysis output, but with less information (Jet_Truth_ID for example). Copies over only basice kinematic information required for the JES study. Was used to simply hadd multiple NTuples during a transition in the analysis procedure.