#!/bin/bash
for beta in 0.440687
do
	for L in 4 8 16 32 64 128 256
	do
		for ini in 1
		do
			python isingmodel.py  $ini $beta $L
		done
	done
done