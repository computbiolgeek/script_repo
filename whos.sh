#!/bin/bash

if [ $# -eq 0 ]; then
  echo "Usage: whos id ..." 1>&2
  exit 1
fi

for id; do
  awk -F : '{print $1, $5}' /etc/passwd | grep -i $id
done
