#!/bin/bash

for value in {1..10}
do 
	#n=$((10000*$value))
	n=100000
	l=$value
	./IE525Project $n $l
done