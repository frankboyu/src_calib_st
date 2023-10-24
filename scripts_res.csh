#!/bin/bash

while read run;
do
  echo $run
  root -b -l -q Resolution.C\(${run}\)
done < runs.dat
