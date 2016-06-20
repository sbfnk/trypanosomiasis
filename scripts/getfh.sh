#!/bin/sh

if [ -z "$1" ]; then exit; fi

rm out/b_$1
touch out/b_$1
for I in ddr ddr_red dd dd_red ddh ddh_red
do
  ./get_b.sh $I""_$1 b_$I >> out/b_$1
done

cat out/b_$1 fullhuman.R | R --no-save
