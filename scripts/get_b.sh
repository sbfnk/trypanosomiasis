#!/bin/sh
pref=""
if [ -n "$2" ];
then
    pref="$2 <- "
fi
    
echo $pref""c\($(tail -n 1 $1 | sed 's/^.* b= //g' | sed 's/ /,/g' | sed 's/[ ,]*$//g')\)
