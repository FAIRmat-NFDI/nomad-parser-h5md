#!/bin/sh

if ! command -v rsync &> /dev/null; then
  echo "rsync required, but not installed!"
  exit 1
else
  rsync -avh nomad-parser-h5md/ .
  rm -rfv nomad-parser-h5md
fi
