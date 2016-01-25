path=/afs/cern.ch/work/n/nmirman/private/CMSSW_7_4_7/src/metsig/Ntuples/Zmumu
date="20160124"

mkdir $path/$date

mkdir $path/$date/Data
mkdir $path/$date/DYJetsToLL
mkdir $path/$date/TTJets
mkdir $path/$date/WJetsToLNu
mkdir $path/$date/WW
mkdir $path/$date/WZ
mkdir $path/$date/ZZ

#cp -i crab_20151025_v1/Run2015D-05Oct2015-v1/res/ntuple_* $path/$date/Data
#cp -i crab_20151025_v1/Run2015D-PromptReco-v4/res/ntuple_* $path/$date/Data
cp crab_20160124/Run2015D-PromptReco-v3/res/ntuple_* $path/$date/Data
cp crab_20160124/DYJetsToLL/res/ntuple_* $path/$date/DYJetsToLL
cp crab_20160124/TTJets/res/ntuple_* $path/$date/TTJets
cp crab_20160124/WJetsToLNu/res/ntuple_* $path/$date/WJetsToLNu
#cp crab_20151025_v1/WW/res/ntuple_* $path/$date/WW
#cp crab_20151025_v1/WZ/res/ntuple_* $path/$date/WZ
#cp crab_20151025_v1/ZZ/res/ntuple_* $path/$date/ZZ
