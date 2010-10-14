#!/bin/sh
echo c\($(tail -n 1 $1 | sed 's/^.* b= //g' | sed 's/ /,/g' | sed 's/[ ,]*$//g')\)
