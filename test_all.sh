#!/bin/sh

param=""

for xi in 0 1 5 10 100 1000
do
    bin/reservoirs $param -l 100 -n 10000 --G._palpalis_gambiense-xi $xi -c $xi,random,none -C xi,groups,habitat 
    param="-d"
    bin/reservoirs -d -l 100 -n 10000 --G._palpalis_gambiense-xi $xi --groups "0;1;2;3;4;5;6;7;8;9;10;11" -c $xi,single,none
    bin/reservoirs -d -l 100 -n 10000 --G._palpalis_gambiense-xi $xi --groups "0;1,2,3;4,5,6,7,8,9,10,11" -c $xi,hum_dom_wild,none
    bin/reservoirs -d -l 100 -n 10000 --G._palpalis_gambiense-xi $xi --groups "0,1,2,3;4,5,6,7,8,9,10,11" -c $xi,humdom_wild,none
    bin/reservoirs -d -l 100 -n 10000 --G._palpalis_gambiense-xi $xi --groups "0;1,2,3,4,5,6,7,8,9,10,11" -c $xi,hum_domwild,none
    bin/reservoirs -d -l 100 -n 10000 --G._palpalis_gambiense-xi $xi -a b -c $xi,random,binary
    bin/reservoirs -d -l 100 -n 10000 --G._palpalis_gambiense-xi $xi -a b --groups "0;1;2;3;4;5;6;7;8;9;10;11" -c $xi,single,binary
    bin/reservoirs -d -l 100 -n 10000 --G._palpalis_gambiense-xi $xi -a b --groups "0;1,2,3;4,5,6,7,8,9,10,11" -c $xi,hum_dom_wild,binary
    bin/reservoirs -d -l 100 -n 10000 --G._palpalis_gambiense-xi $xi -a b --groups "0,1,2,3;4,5,6,7,8,9,10,11" -c $xi,humdom_wild,binary
    bin/reservoirs -d -l 100 -n 10000 --G._palpalis_gambiense-xi $xi -a b --groups "0;1,2,3,4,5,6,7,8,9,10,11" -c $xi,hum_domwild,binary
    bin/reservoirs -d -l 100 -n 10000 --G._palpalis_gambiense-xi $xi -a f -c $xi,random,fractional
    bin/reservoirs -d -l 100 -n 10000 --G._palpalis_gambiense-xi $xi -a f --groups "0;1;2;3;4;5;6;7;8;9;10;11" -c $xi,single,fractional
    bin/reservoirs -d -l 100 -n 10000 --G._palpalis_gambiense-xi $xi -a f --groups "0;1,2,3;4,5,6,7,8,9,10,11" -c $xi,hum_dom_wild,fractional
    bin/reservoirs -d -l 100 -n 10000 --G._palpalis_gambiense-xi $xi -a f --groups "0,1,2,3;4,5,6,7,8,9,10,11" -c $xi,humdom_wild,fractional
    bin/reservoirs -d -l 100 -n 10000 --G._palpalis_gambiense-xi $xi -a f --groups "0;1,2,3,4,5,6,7,8,9,10,11" -c $xi,hum_domwild,fractional
done
