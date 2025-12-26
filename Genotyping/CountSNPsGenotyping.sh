#!/usr/bin/env bash

#Counting SNPs

for i in {1..2}; do cut -c $i Calls.Comp.geno | grep -v 9 | wc -l; done