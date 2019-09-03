#!/bin/bash
outfile=RooUnfoldFitExample.py.bbb.ref
python -u examples/RooUnfoldFitExample.py bbb --mode RooUnfoldSpec > $outfile
bash ref/cleanup.sh $outfile
diff $outfile ref/$outfile
