path=/afs/cern.ch/work/n/nmirman/private/CMSSW_7_4_7/src/metsig/Ntuples/Zmumu
date="20160418"

mkdir $path/$date

mkdir $path/$date/Data
mkdir $path/$date/DYJetsToLL
mkdir $path/$date/TTJets
mkdir $path/$date/WJetsToLNu
mkdir $path/$date/WW
mkdir $path/$date/WZ
mkdir $path/$date/ZZ
mkdir $path/$date/ST_tW_top
mkdir $path/$date/ST_tW_antitop

cp crab_$date/Run2015D/res/ntuple_* $path/$date/Data
cp crab_$date/DYJetsToLL/res/ntuple_* $path/$date/DYJetsToLL
cp crab_$date/TTJets/res/ntuple_* $path/$date/TTJets
cp crab_$date/WJetsToLNu/res/ntuple_* $path/$date/WJetsToLNu
cp crab_$date/WW/res/ntuple_* $path/$date/WW
cp crab_$date/WZ/res/ntuple_* $path/$date/WZ
cp crab_$date/ZZ/res/ntuple_* $path/$date/ZZ
cp crab_$date/ST_tW_top/res/ntuple_* $path/$date/ST_tW_top
cp crab_$date/ST_tW_antitop/res/ntuple_* $path/$date/ST_tW_antitop
