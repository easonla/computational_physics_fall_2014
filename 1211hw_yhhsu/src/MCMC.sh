#!/bin/bash
for beta in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
do
	for L in 128 256
	do
		for ini in 0 1
		do
			python isingmodel.py  $ini $beta $L
		done
	done
done
