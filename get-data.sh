#!/bin/bash

#common location on server
loc='/Volumes/data/gac_suncor/hsay/gac_samples/allvall_plasmids'

# edit these for each sample
dir='2-T4_may14-21'
dir_lo='2-T4'

# make a directory for the new data
mkdir /Users/pemadorjee/Documents/thesis/data/$dir_lo

# then get the data
scp pdorjee@129.100.236.35:$loc/$dir/initial_polished/bakta/*.tsv data/$dir_lo/
