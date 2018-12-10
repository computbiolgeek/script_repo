#!/bin/bash

if [ $# -ne 2 ]; then
  echo "Usage: $0 file1 file2" 1>&2
  echo "file1: list of correct spellings" 1>&2
  echo "file2: file to be checked" 1>&2
fi

if [ ! -r $1 ]; then
  echo "$0: $1 is not readable" 1>&2
  exit 1
fi

if [ ! -r $2 ]; then
  echo "$0: $2 is not readable" 1>&2
  exit 1
fi

aspell list < $2 |
while read line; do
  if ! grep -qw "$line" $1 > /dev/null; then
    echo $line
  fi
done
