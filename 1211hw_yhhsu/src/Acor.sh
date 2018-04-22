#!/bin/bash
for beta in 0.1 0.2 0.3 0.4 0.5 1.0
do
	for L in 4 8 16 32 64
	do
		for ini in 0
		do
			python acor.py  $ini $beta $L
		done
	done
done
