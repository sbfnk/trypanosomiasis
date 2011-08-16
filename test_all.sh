#!/bin/sh

bin/reservoirs -l 100 -n 100000 > single_none.csv
bin/reservoirs -l 100 -n 100000 --groups "0;1,2,3;4,5,6,7,8,9,10,11" > hum_dom_wild_none.csv
bin/reservoirs -l 100 -n 100000 --groups "0,1,2,3;4,5,6,7,8,9,10,11" > humdom_wild_none.csv
bin/reservoirs -l 100 -n 100000 --groups "0;1,2,3,4,5,6,7,8,9,10,11" > hum_domwild_none.csv
bin/reservoirs -l 100 -n 100000 -a b > single_binary.csv
bin/reservoirs -l 100 -n 100000 -a b --groups "0;1,2,3;4,5,6,7,8,9,10,11" > hum_dom_wild_binary.csv
bin/reservoirs -l 100 -n 100000 -a b --groups "0,1,2,3;4,5,6,7,8,9,10,11" > humdom_wild_binary.csv
bin/reservoirs -l 100 -n 100000 -a b --groups "0;1,2,3,4,5,6,7,8,9,10,11" > hum_domwild_binary.csv
bin/reservoirs -l 100 -n 100000 -a f > single_fractional.csv
bin/reservoirs -l 100 -n 100000 -a f --groups "0;1,2,3;4,5,6,7,8,9,10,11" > hum_dom_wild_fractional.csv
bin/reservoirs -l 100 -n 100000 -a f --groups "0,1,2,3;4,5,6,7,8,9,10,11" > humdom_wild_fractional.csv
bin/reservoirs -l 100 -n 100000 -a f --groups "0;1,2,3,4,5,6,7,8,9,10,11" > hum_domwild_fractional.csv
