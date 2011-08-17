#!/bin/sh

bin/reservoirs -l 100 -n 100000 -C groups,habitat -c single,none
bin/reservoirs -l 100 -n 100000 --groups "0;1,2,3;4,5,6,7,8,9,10,11" -c hum_dom_wild,none
bin/reservoirs -l 100 -n 100000 --groups "0,1,2,3;4,5,6,7,8,9,10,11" -c humdom_wild,none
bin/reservoirs -l 100 -n 100000 --groups "0;1,2,3,4,5,6,7,8,9,10,11" -c hum_domwild,none
bin/reservoirs -l 100 -n 100000 -a b -c single,binary
bin/reservoirs -l 100 -n 100000 -a b --groups "0;1,2,3;4,5,6,7,8,9,10,11" -c hum_dom_wild,binary
bin/reservoirs -l 100 -n 100000 -a b --groups "0,1,2,3;4,5,6,7,8,9,10,11" -c humdom_wild,binary
bin/reservoirs -l 100 -n 100000 -a b --groups "0;1,2,3,4,5,6,7,8,9,10,11" -c hum_domwild,binary
bin/reservoirs -l 100 -n 100000 -a f -c single,fractional
bin/reservoirs -l 100 -n 100000 -a f --groups "0;1,2,3;4,5,6,7,8,9,10,11" -c hum_dom_wild,fractional
bin/reservoirs -l 100 -n 100000 -a f --groups "0,1,2,3;4,5,6,7,8,9,10,11" -c humdom_wild,fractional
bin/reservoirs -l 100 -n 100000 -a f --groups "0;1,2,3,4,5,6,7,8,9,10,11" -c hum_domwild,fractional
