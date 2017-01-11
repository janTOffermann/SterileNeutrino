# SterileNeutrino

### Summary
This code, consisting of the script DetectorSim3N3.C, was created to simulate a proposed experiment for detecting keV-mass sterile neutrinos. Using information from the project proposal as well as the Chart of Nuclides (from Brookhaven National Laboratory), the script simulates K-capture events for Cesium-131 atoms, and simulates the produced Xenon-131 recoil ions and Auger electrons as they propagate into multi-channel plate detectors, passing through constant electric and magnetic fields. The simulation stores information from the detectors, and provides a function to reconstruct a histogram of the square-mass of the neutrinos produced.

### Usage
The script can be broken down into two main functions. One of the functions, "Sim", runs the simulation while a second function, "GenHist", can be used to analyze the results as they are saved in a .root file. Below are the syntaxes for these two functions.
```sh
void Sim(UInt limit, Double_t mixAngle, Bool_t no_l_shell = kTRUE, Bool_t ERRORS = kFALSE)
void GenHist(TString fname_root = "sim.root")
```
  - limit is the number of K-capture events to run for (unsigned int)
  - mixAngle is the % of events to which to assign massive neutrinos (double)
  - no_l_shell determines whether or not to simulate events with L-shell gamma rays (bool, default is true)
  - ERRORS determines whether or not to include errors resulting from detector precisions (bool, default is false)
  - fname_root is the name of the .root file containing the results from Sim (string, default is "sim.root")

It is important to note that the massive neutrinos have a hard-coded mass value of 10 keV - this is a global variable which would need to be changed in the code itself.
As errors in recoil ion position appeared to have drastic results (the square-mass histograms produced with these errors showed no clear peaks and were not valuable for analysis), the lines of code that add errors to recoil ion position have been commented out by default.

### Dependencies
[![N|Solid](https://d35c7d8c.web.cern.ch/sites/d35c7d8c.web.cern.ch/files/website-banner-allnew-croped_3.png)](https://root.cern.ch)

This script makes use of ROOT, the data analysis framework developed at CERN (see https://root.cern.ch for information on installing ROOT).
The code was run in the ROOT interpreter. For example, one may run the following (if renaming the .root file from the default of "sim.root" after running Sim):
```sh
$ root -l
$ root [0] .L DetectorSim3N3.C+
$ root [1] Sim(50000, 0.0001, kTRUE, kTRUE)
$ root [2] GenHist("50k_events_with_errors.root")
```
### Acknowledgements

I would like to thank Professor Peter Meyers for serving as my advisor for this project. I would also like to thank my father Edmond Offermann for helping me make the most of ROOT's functionality, and helping me hunt down various bugs during the development of this code.

### Contact
For any questions, I may be reached at jano@princeton.edu. Thank you.

Jan Tuzlic Offermann
