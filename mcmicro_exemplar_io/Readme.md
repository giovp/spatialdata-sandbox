To download the data and generate the output of the pipeline do the following:

```
mkdir data
cd data

# download the data
nextflow run labsyspharm/mcmicro/exemplar.nf --name exemplar-001 --path .
nextflow run labsyspharm/mcmicro/exemplar.nf --name exemplar-002 --path .

# run mcmicro
# the first command takes 5-10 min, the second takes 30-40 min
# see more here: https://mcmicro.org/tutorial/tutorial.html
nextflow run labsyspharm/mcmicro --in exemplar-001
nextflow run labsyspharm/mcmicro --in exemplar-002
```
