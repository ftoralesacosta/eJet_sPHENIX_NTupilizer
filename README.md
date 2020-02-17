
# Getting Started
This code is adapted from the JES scale code, which is a fork of the sPHENIX/macros repository and edits Jin Huang's "myjetanalysis" code.

This code generates the Jet Energy Scale (JES) correction factor, as well as parametrizes the Jet Energy Resolution (JER).
The code that does this is in the JES directory. The code needs to be compiled, but is compiled separatly from Fun4All, and simply needs ROOT installed. All that needs to be done is to point it at the correct NTuple. The following is only required is you want make your own MonteCarlo production.

The code used for fun4all is adapted from https://github.com/sPHENIX-Collaboration/macros, and works on TTrees made from https://github.com/sPHENIX-Collaboration/tutorials/tree/master/myjetanalysis 

The myjetanalysis code needs to built and installed with Fun4All, as outlined here: https://wiki.bnl.gov/sPHENIX/index.php/Example_of_using_DST_nodes

The code was originally run inside a Singularity container, which can be found here:
https://github.com/sPHENIX-Collaboration/Singularity

# Running

Once Fun4All with the myjetanalysis is built, and the JES code is compiled with the included Makefile, simply run by:

`./JES_myjetanalysis [myjetanalysis_output.root]`

There are two main loops to the JES analysis. The first loop goes through reconstructed and truth jets, and fits gausians to the distributions in pT Truth bins to eventualy obtain the correction as a function of Reconstructed pT (saved as a TF1). The second loop then goes through reconstructed jets to apply the correction.

The result will be an output root file with several histograms/TGraphs, and the TF1 correction function. Several PDFs of plots will also be generated. Apologies for any clutter.

All the machinery is in place to run your own Monte Carlo Production. The main Fun4All_G4_sPHENIX.C macro calls the myjetanalysis code to create TTrees with jet variables.

# Changing Bins and Adding Dependencies

The bins are defined with arrays starting on Line 625. Simply change the min and max Truth pT, and the bin width (currently 2.5 GeV/c).

By default, the only binning is in Jet pT, however adding eta dependency, or z-vertex dependance should be relativley simple. The implementation uses unfortunatly clunky strings to ensure unique historgram names, but is working as intended. I wrote eta depenence as an example in Line 651. The steps to add additional dependencies are as the follows:

1. Create a vector for the the bin edges and a vector for the bin centers (shown in Line 651)
2. Add an extra dimension to N-dimensional TH1F vector `v_Eta_pT_TH1F` and the N-dimensional TH2F vector `v_Eta_pT_TH2F`
3. Add an extra dimension to the string vectors `root_cut_string` and `title_cut_string`
4. Fill the string vectors shown in Line 672-673, using the newly created vectors containing bin information.
5. Finally, add additional for loops for the new bins for the reconstructed and correction loops (shown in Lines 728 & 840. All the functions inside the loops take a string in 1-D TH1F vector as arguments. The outer loop is simply changed to go throught the N-dimensional vectors. In this way none of the functions need editing if all one wants to do is add additional dependencies.

# To Do:
Hard coded binning is of course a pain. Adding TEnv support or some external configuration file would be great.

*_Large MC prodcution_*. The last iteration using this code had a hard cut on lower Jet pT, and utilized MC based on a much older version of Fun4All. An new large scale production implementing PhaseSpace:bias2SelectionPow in the pythia configuration file is highly recomended. The current pythia configuration file https://github.com/ftoralesacosta/macros/blob/master/macros/g4simulations/phpythia8.cfg has this implemented. The power and reference pT hat could use some tweaking.
