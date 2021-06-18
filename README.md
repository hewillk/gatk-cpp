## Introduction
Hardware-friendly version of GATK4 HaplotypeCaller.

## Getting started

    mkdir build
    cd build
    cmake ..
    make
    ./gatk -I ../input/chrM.sam -R ../reference/chrM.fa -O ../output/chrM.vcf 
