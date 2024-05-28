# cfp-optimizer

Carbon footprint minimization through workload management

## Optimized job scheduling across multiple data centers with green objectives

The greenness of the power sources to data centers vary by location and time. Real-time data about the power source mix to data centers is becoming increasingly available, and future predictions are also viable. This provides an opportunity to minimize the carbon footprint due to computing workloads by properly scheduling jobs through a choice of a data center and a particular time of day.

We consider a cloud computing infrastructure consisting of multiple data centers that is subject to a batch job workload. When submitting a job, the user specifies a deadline for the start time of the job. A job declares its resource requirements and an estimate of its run time. Also, jobs may have different priorities. There may be other constraints related to a job executing on particular data centers. A job scheduler makes a decision for each job on (1) an assigned data center and (2) a start time for the job.
The objective of the scheduler is to minimize the overall carbon footprint due to processing the workload, while satisfying the constraints.

We consider two cases: off-line and on-line scheduling. In the off-line case, the optimization problem is solved once for a given set of jobs over a time window in the future. And, in the on-line case, the optimization problem is solved periodically, potentially adjusting prior solutions.
