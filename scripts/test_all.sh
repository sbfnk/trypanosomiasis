#!/bin/sh

bin/reservoirs -l 100 -s 100000 -c random,none -C groups,habitat 
bin/reservoirs -d -l 100 -s 100000 --groups "0;1;2;3;4;5;6;7;8;9;10;11" -c single,none
bin/reservoirs -d -l 100 -s 100000 --groups "0;1,2,3;4,5,6,7,8,9,10,11" -c hum_dom_wild,none
bin/reservoirs -d -l 100 -s 100000 --groups "0,1,2,3;4,5,6,7,8,9,10,11" -c humdom_wild,none
bin/reservoirs -d -l 100 -s 100000 --groups "0;1,2,3,4,5,6,7,8,9,10,11" -c hum_domwild,none
bin/reservoirs -d -l 100 -s 100000 -a b -c random,binary
bin/reservoirs -d -l 100 -s 100000 -a b --groups "0;1;2;3;4;5;6;7;8;9;10;11" -c single,binary
bin/reservoirs -d -l 100 -s 100000 -a b --groups "0;1,2,3;4,5,6,7,8,9,10,11" -c hum_dom_wild,binary
bin/reservoirs -d -l 100 -s 100000 -a b --groups "0,1,2,3;4,5,6,7,8,9,10,11" -c humdom_wild,binary
bin/reservoirs -d -l 100 -s 100000 -a b --groups "0;1,2,3,4,5,6,7,8,9,10,11" -c hum_domwild,binary
bin/reservoirs -d -l 100 -s 100000 -a f -c random,fractional
bin/reservoirs -d -l 100 -s 100000 -a f --groups "0;1;2;3;4;5;6;7;8;9;10;11" -c single,fractional
bin/reservoirs -d -l 100 -s 100000 -a f --groups "0;1,2,3;4,5,6,7,8,9,10,11" -c hum_dom_wild,fractional
bin/reservoirs -d -l 100 -s 100000 -a f --groups "0,1,2,3;4,5,6,7,8,9,10,11" -c humdom_wild,fractional
bin/reservoirs -d -l 100 -s 100000 -a f --groups "0;1,2,3,4,5,6,7,8,9,10,11" -c hum_domwild,fractional
