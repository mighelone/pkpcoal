# Comparison of performances

## Anaconda + Numba

>>> time ./runPKP input.yml -n 2

real	0m26.614s
user	0m47.076s
sys	0m1.158s

## Python + Numpy

>>> time ./runPKP input.yml -n 2

real	0m38.399s
user	1m9.661s
sys	0m1.275s
