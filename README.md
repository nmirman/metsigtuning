# METSigTuning
Optimization code and ntuple maker for updated MET significance algorithm, based on the 13 TeV CMS dataset.  Code is interfaced with the CMS Software (CMSSW).

### MakeNtuples
Runs over CMS datasets to select events for analysis containing at least two muons.  Events are stored in ntuple format that is readable in `MetSigFit`.

### METSigFit
Updated version of MET significance algorithm, including a maximum likelihood fit for parameter optimization.  Computes a covariance matrix corresponding to the missing transverse energy in each event.  Estimates the statistical significance of the MET variable.
