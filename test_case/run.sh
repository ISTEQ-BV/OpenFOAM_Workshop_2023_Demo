#!/bin/sh

set -eu

series="my_machine-4k"

results_dir="benchmarking/$series"
mkdir -p "$results_dir"

#for n in 15 16 17 18 19 20  22  23  25  27  29  32  34  37  40  43  47  50  54  59  63  68  74  80  86  93 100 108 117 126 136 147 159 172 185 200
#for n in 16 17 18 19 20  22  23  25  27  29  32  34  37  40  43  47  50  54  59  63  68  74  80  86  93 100 108 117 126 136
#for n in 16 17 18 19 20  22  23  25  27  29  32  34  37  40  43  47  50  54  59  63  68  74  80  86  93 100 108 117
do
    rm -rf 1 constant/cellDist
    sed "s/NNN/$n/g" system/blockMeshDict.template > system/blockMeshDict
    blockMesh
    renumberMesh -dict system/renumberMeshDict -constant
    field_traversal_benchmark
    mv results.csv "$results_dir/results-$n.csv"
done

