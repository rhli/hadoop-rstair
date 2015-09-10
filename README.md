Overview
=====

In clustered file systems, single node failure are most common, however, there
are failure bursts happen from time to time.  We thus propose RSTAIR code for a
fast recovery for single node failure, while tolerable of burst of failures
with low storage overhead.

User's Guide
=====
R-STAIR code is implemented and integrated with Facebook's Hadoop-20.

I. Install Hadoop-20
====

Hadoop-20 can be downloaded from [here](https://github.com/facebookarchive/hadoop-20).  
Download and extract it in one of your hosts.
If the above link is not working, you can download from [homepage of this project] (http://ansrlab.cse.cuhk.edu.hk/software/ear).

To simplify our description, in this instruction, we refer to the hadoop-20/ directory as hadoop directory.  

Call the following command to update the environment variable.
> export HADOOP_HOME=*absolute address of hadoop directory*   
> export PATH=$PATH:$HADOOP_HOME/bin   

It is a good practice to include the above two lines in your ~/.bashrc file.

Add the following line to $HADOOP_HOME/conf/hadoop-env.sh
> export JAVA_HOME=*YOUR JAVA_HOME*

II. Installation and Configuration
====

0. Run `bash install.sh` to install R-STAIR code.

1. Configure masters and slaves  
  To start with, you need to make sure in your cluster every machine is
  public-key accessible from other machines.  
  Input the IP address/hostname of master and slaves in file conf/masters and
  conf/slaves.  

2. Create cluster topology definition script.
  Hadoop uses a topology definition script for rack awareness.  
  The script takes IP address or host name of a machine as input and returns
  which rack the machine belongs to.  Take a look at $HADOOP_HOME/rackAware.sh
  to get a feeling and change it according to your cluster topology.

  To use R-STAIR code, we assume you have enough racks and enough nodes in each
  rack.  For example, if you set n=6 and r=5, you should have at least 6 racks
  and at least 5 nodes per rack.

3. Configure conf/core-site.xml conf/hdfs-site.xml conf/raid.xml   
  according to the instructions in the xml files.

4. Finishing: 
  Copy the hadoop directory to every machine in the same location.

III. Configuration
====
To change parameters of R-STAIR code, please refer to 
$HADOOP_HOME/src/contrib/raid/src/java/raid-default.xml as an example.

To make R-STAIR code work properly, please also update corresponding
configurations in $HADOOP_HOME/conf/hdfs-site.xml for the R-STAIR code
placement.

Developer's Guide
=====

Our modifications can be divided into three parts:

1. Coding implementation 
    R-STAIR is implemented in C++ and integrated with HDFS-RAID through JNI.
    Find our implementation in src/native/src/org/apache/hadoop/util/stair.c
2. R-STAIR code placement
    We implement our two-dimension round-robin placement to place data/parity
    blocks to achieve load balance and exploit the locality property of R-STAIR
    codes.  (src/contrib/raid/src/java/org/apache/hadoop/hdfs/server/namenode)
      - See BlockPlacementPolicyStair.java for more detail.
3. MapReduce integration
    We explore how to leverage the feature of R-STAIR code to improve MapReduce
    performance.  We modify MapReduce scheduler to integrate 
      - We mainly modified JobInProgress.java to leverage locality of R-STAIR
        code for MapReduce jobs.

Change Log
=====
Version 1.0.0 (September, 2015): First release. 

