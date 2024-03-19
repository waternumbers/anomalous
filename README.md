# anomalous

R implimentations of the pelt and capa algorithms for changepoint and anomaly detection.

Implimnetation is build around an S3 class for partitions and R6 classes for cost functions. 
The capa and pelt implimentatons make use of the methods of the partition class. This allows for the replacement of the current partition class with a new class with similarly names methods.

The partition class makes us of the cost classes via names methods, again allowing for generalisation.

## TO DO

- Handle missing values in the following cost functions
  - local_reg 
  - lad
  - gauss_reg
