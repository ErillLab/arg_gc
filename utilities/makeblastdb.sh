#!usr/bin/bash

makeblastdb -in $1 -out $2 -dbtype 'prot'

rm $1
