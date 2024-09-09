** THIS IS SOFTWARE IN DEVELOPMENT - EXPECT BUGS AND CHANGES **

# anomalous

R implementations of the pelt and capa algorithms for changepoint and anomaly detection.

The implementation is build around an S3 class for partitions and R6 classes for cost functions. 
The capa and pelt implementations make use of the methods of the partition class. This allows for the replacement of the current partition class with a new class with similarly names methods.

The partition class makes us of the cost classes via names methods, again allowing for generalisation.


## TO DO

- Documentation
- Check handling of missing data
- Consistency in the input formats to the cost function
- Improve plotting methods especially for multivariate data
- Add handling of multivariate anomaly detection
