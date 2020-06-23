# Platypus germline pipeline
Author: Nagarajan Paramasivam, DKFZ

Pipeline calls variants from single or multiple WES/WGS from an family and annotates the variants with gene and allele frequency information.
# 
# To Run the pipeline,

1. For single sample, the results are written to a subdir respect to the PID

```
# set CANCERorGERMLINE=2 in config file
sh run_platypusPipeline.sh -c CONFIG_FILE $PID_1
sh run_platypusPipeline.sh -c CONFIG_FILE $PID_2
```
2. This will do a multi calling, the results are written to a PID_1 subdir.
```
#set CANCERorGERMLINE=2
sh run_platypusPipeline.sh -c CONFIG_FILE $PID_1 $PID_2 $PID_3 ... $PID_Inf
#note: if there are Inf samples, it is better to go for single calls, and do a joint genotyping. 
```
3. For cancer, for blood and tumor
```
#set CANCERorGERMLINE=1
sh run_platypusPipeline.sh -c CONFIG_FILE $PID_1
```
