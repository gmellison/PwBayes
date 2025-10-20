#! /bin/bash

for ii in 1 2 3
  do
  nohup Rscript sim2_altfossil.R $ii & # run in background
  done

