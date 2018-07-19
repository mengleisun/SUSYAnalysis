ls ../../DATA/SingleMuon/|xargs -I{} -P 4 root -b -q "analysis_hadron.C+(\"{}\")"
