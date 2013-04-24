#!/bin/bash

for i in 148 149 150 151 170 171 172 173 178 179 180 181 182 183 184 185 204 237
do
	mpirun -n 8 ./run$i/HIM < ./run$i/temp_in
done






