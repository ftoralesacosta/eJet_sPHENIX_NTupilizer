# Extremely Interesting Background Information

This code is adapted from the JES scale code, which I orginally forked from the sPHENIX/macros repository and edits Jin Huang's "myjetanalysis" code. The code now relies on a separate cloning the the sPHENIX/macros repo, but is still maintained as a fork of tutorials/myjetanalysis code.

This code simulates electron proton collisions using the sPHENIX detector, and uses hybrid tracks. It outputs simple NTuples with true electron information, as well as reconstructed and truth level jets. Rey Cruz Torres has implemented an **all silicon** tracking design for a future EIC based on the same code based (denoted with AllSi prefixes) that outputs NTupples with same structure.

(Depricated as of 7/2020): The code in the JES directory generates the Jet Energy Scale (JES) correction factor, as well as parametrizes the Jet Energy Resolution (JER). The code needs to be compiled, but is compiled separatly from Fun4All, and simply needs ROOT installed. All that needs to be done is to point it at the correct NTuple if simulations are provided. 

# Getting Started
1. The code was originally ran inside a Singularity container, which can be found here:
https://github.com/sPHENIX-Collaboration/Singularity
One needs to run ./updatebuild.sh and follow the steps in the README of that repo.

2. One needs to clone the sPHENIX macro repository: https://github.com/sPHENIX-Collaboration/macros if they want to run their own simulation.

3. The myjetanalysis code needs to built and installed with Fun4All, as outlined here: https://wiki.bnl.gov/sPHENIX/index.php/Example_of_using_DST_nodes. References to <sourcedir> for this repository mean `e_Jet_sPHENIX/myjetanalysis/src/ `, where one should see the autogen.sh file. Create a `build` and `install` directory (I suggest in the same directory that holds this repo) and follow the instructions under the "Building a package" section from the link.
  
4. Compile the JES code separately. There is a JES directory with an included make file. Once compiled, just run with:
`./e_Jet_Analysis.cc <Input NTupple>` for simple electron+jet kinematic analysis. The JES is in progress from pp and AA collisions settings.

# Event Production
First, [compile the code](https://wiki.bnl.gov/sPHENIX/index.php/Example_of_using_DST_nodes) and then make sure to source `local_setup.sh $MYINST` so the compiled libraries can be used in the macros

## Copy the PYTHIA8 Configuration for electron-proton collisions
Copy the phpythia8.cfg file into the macros/macros/g4simulation/ directory (note: this will overwrite the default file)

## The Fun4All_G4_sPHENIX.C and G4_Jets.C macro must be edited. 
They can be found in macros/macros/g4simulation/ in the sPHENIX macros repository (https://github.com/sPHENIX-Collaboration/macros).

### Fun4All_G4_sPHENIX.C

Line 31 (right before the first R_LOAD_LIBRARY):
  ```
  #include <myjetanalysis/MyJetAnalysis.h>
  R__LOAD_LIBRARY(libmyjetanalysis.so)
  ```
Line 70:
     const bool runpythia8 = false; -> true

Line 81:
     const bool particles = true && !readhits; -> false

Line 601 (right before event processing):
  ```
  gSystem->Load("libmyjetanalysis");
  std::string jetoutputFile = std::string(outputFile) + std::string("electrons+jets.root");
  MyJetAnalysis *myJetAnalysis = new MyJetAnalysis("AntiKt_Track_r10","AntiKt_Truth_r10",jetoutputFile.data());
  se->registerSubsystem(myJetAnalysis);
  ```
### G4_Jets.C

Line 33:
 ```
 truthjetreco->add_algo(new FastJetAlgo(Jet::ANTIKT,1.0),"AntiKt_Truth_r10");
 ```
Line 81:
```
trackjetreco->add_algo(new FastJetAlgo(Jet::ANTIKT,1.0),"AntiKt_Track_r10");

```

## Reco-Truth Jet Matching
The original source code uses the function *max_truth_jet_by_enengy(jet)*.
This has been changed to *unique_truth_jet_from_reco(jet);* in order to obtain 1:1 truth-reco matches. This is done at the exclusion of fake  jets in order to specifically study the jet reconstruction efficency, Energy/Momentum scale, and angular resolutions where a precise matching is important. For future studies of the fake rate of jets, one should likely switch to the original matching function.

# Running Jet Energy Resolution

Once Fun4All with the myjetanalysis is built, and the JES code is compiled with the included Makefile, simply run by:

`./JES_myjetanalysis [myjetanalysis_output.root]`

There are two main loops to the JES analysis. The first loop goes through reconstructed and truth jets, and fits gausians to the distributions in pT Truth bins to eventualy obtain the correction as a function of Reconstructed pT (saved as a TF1). The second loop then goes through reconstructed jets to apply the correction.

The result will be an output root file with several histograms/TGraphs, and the TF1 correction function. Several PDFs of plots will also be generated. Apologies for any clutter.

All the machinery is in place to run your own Monte Carlo Production. The main Fun4All_G4_sPHENIX.C macro calls the myjetanalysis code to create TTrees with jet variables.

<!---
# Changing Bins and Adding Dependencies

The bins are defined with arrays starting on Line 625. Simply change the min and max Truth pT, and the bin width (currently 2.5 GeV/c).

By default, the only binning is in Jet pT, however adding eta dependency, or z-vertex dependance should be relativley simple. The implementation uses unfortunatly clunky strings to ensure unique historgram names, but is working as intended. I wrote eta depenence as an example in Line 651. The steps to add additional dependencies are as the follows:

1. Create a vector for the the bin edges and a vector for the bin centers (shown in Line 651)
2. Add an extra dimension to N-dimensional TH1F vector `v_Eta_pT_TH1F` and the N-dimensional TH2F vector `v_Eta_pT_TH2F`
3. Add an extra dimension to the string vectors `root_cut_string` and `title_cut_string`
4. Fill the string vectors shown in Line 672-673, using the newly created vectors containing bin information.
5. Finally, add additional for loops for the new bins for the reconstructed and correction loops (shown in Lines 728 & 840. All the functions inside the loops take a string in 1-D TH1F vector as arguments. The outer loop is simply changed to go throught the N-dimensional vectors. In this way none of the functions need editing if all one wants to do is add additional dependencies.
-->

# To Do:
Hard coded binning is of course a pain. Adding TEnv support or some external configuration file would be great.

*_Large MC prodcution_*. The last iteration using this code had a hard cut on lower Jet pT, and utilized MC based on a much older version of Fun4All. An new large scale production implementing PhaseSpace:bias2SelectionPow in the pythia configuration file is highly recomended. The current pythia configuration file https://github.com/ftoralesacosta/macros/blob/master/macros/g4simulations/phpythia8.cfg has this implemented. The power and reference pT hat could use some tweaking.
