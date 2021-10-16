#!/bin/bash

# Seems silly to integrate this with check_efit_results.sh since it is only used in this test
# but that could be done...

if ! test -e fitout.dat; then
  echo "File does not exist: fitout.dat"
  exit 1
fi
