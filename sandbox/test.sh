#!/bin/sh
if [[ "$1" =~ ^(darter|yellowstone)$ ]]; then
echo haha
else
echo "say wha"
echo $(ls *)
fi
