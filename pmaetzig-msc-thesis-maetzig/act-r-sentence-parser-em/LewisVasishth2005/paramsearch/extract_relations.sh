#!/bin/bash

simset=( `awk < paramsearch.txt '{print $2}' `)

for set in ${simset[@]}
do
    awk < $set-trace.txt '/SENTENCE:|"subject"|"antecedent"/' > $set-relations.txt
done
