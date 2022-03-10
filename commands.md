#!/usr/bin/env python3

# -_- coding: utf-8 -_-

"""
Created on Wed Mar 9 14:38:24 2022

@author: François Clapuyt
"""

# Test Otemma

in_dtm="/home/elic/fclapuyt/SedInConnect_2.3/inputs/Otemma/Otemma_filledWL2019.tif"
in_cs=2
out_w="/home/elic/fclapuyt/SedInConnect_2.3/outputs/Otemma/otemma_weight.tif"
out_ic="/home/elic/fclapuyt/SedInConnect_2.3/outputs/Otemma/otemma_ic.tif"

CavalliConnectivityout(in_dtm, in_cs, in_dtm, out_ic)

```bash
 [Mar/09] fclapuyt@pelican 15:29:33 ~/SedInConnect_2.3 22 > python -i sedinconnect_cli.py
>>> in_dtm="/home/elic/fclapuyt/SedInConnect_2.3/inputs/Otemma/Otemma_filledWL2019.tif"
>>> in_cs=2
>>> out_w="/home/elic/fclapuyt/SedInConnect_2.3/outputs/Otemma/otemma_weight.tif"
>>> out_ic="/home/elic/fclapuyt/SedInConnect_2.3/outputs/Otemma/otemma_ic.tif"
>>> CavalliConnectivityout(in_dtm, in_cs, in_dtm, out_ic)
09/03/2022 15:30:02     opening DTM for Connectivity was successful!
D8FlowDir version 5.3.9
Input file /home/elic/fclapuyt/SedInConnect_2.3/inputs/Otemma/Otemma_filledWL2019.tif has projected coordinate system.
Nodata value input to create partition from file: -99999.000000
Nodata value recast to float used in partition raster: -99999.000000
This run may take on the order of 11 minutes to complete.
This estimate is very approximate.
Run time is highly uncertain as it depends on the complexity of the input data
and speed and memory of the computer. This estimate is based on our testing on
a dual quad core Dell Xeon E5405 2.0GHz PC with 16GB RAM.
All slopes evaluated. 0 flats to resolve.
Processors: 8
Header read time: 1.159405
Data read time: 1.839811
Compute Slope time: 20.121960
Write Slope time: 30.485983
Resolve Flat time: 24.139471
Write Flat time: 10.626566
Total time: 88.373197
09/03/2022 15:31:43     opening fd8 was successful!
DinfFlowDir version 5.3.9
Input file /home/elic/fclapuyt/SedInConnect_2.3/inputs/Otemma/Otemma_filledWL2019.tif has projected coordinate system.
Nodata value input to create partition from file: -99999.000000
Nodata value recast to float used in partition raster: -99999.000000
This run may take on the order of 11 minutes to complete.
This estimate is very approximate.
Run time is highly uncertain as it depends on the complexity of the input data
and speed and memory of the computer. This estimate is based on our testing on
a dual quad core Dell Xeon E5405 2.0GHz PC with 16GB RAM.
All slopes evaluated. 0 flats to resolve.
Processors: 8
Header read time: 1.287074
Data read time: 1.909564
Compute Slope time: 50.665171
Write Slope time: 43.562313
Resolve Flat time: 24.635488
Write Flat time: 34.689705
Total time: 156.749316
09/03/2022 15:34:48     opening sd8 was successful!
09/03/2022 15:34:54     opening Weight was successful!
AreaD8 version 5.3.9
Input file /home/elic/fclapuyt/SedInConnect_2.3/inputs/Otemma/Otemma_filledWL2019p.tif has projected coordinate system.
This run may take on the order of 1 minutes to complete.
This estimate is very approximate.
Run time is highly uncertain as it depends on the complexity of the input data
and speed and memory of the computer. This estimate is based on our testing on
a dual quad core Dell Xeon E5405 2.0GHz PC with 16GB RAM.
Nodata value input to create partition from file: -32768.000000
Nodata value recast to int16_t used in partition raster: -32768
Number of Processes: 8
Read time: 1.902497
Compute time: 45.694717
Write time: 7.605200
Total time: 55.202415
09/03/2022 15:37:05     Process concluded succesfully
 42.67 seconds for Weighted-Flow Length calculation with  9115  iterations
AreaDinf version 5.3.9
Input file /home/elic/fclapuyt/SedInConnect_2.3/inputs/Otemma/Otemma_filledWL2019ang.tif has projected coordinate system.
This run may take on the order of 1 minutes to complete.
This estimate is very approximate.
Run time is highly uncertain as it depends on the complexity of the input data
and speed and memory of the computer. This estimate is based on our testing on
a dual quad core Dell Xeon E5405 2.0GHz PC with 16GB RAM.
Nodata value input to create partition from file: -340282346638528859811704183484516925440.000000
Nodata value recast to float used in partition raster: -340282346638528859811704183484516925440.000000
Processors: 8
Read time: 3.596958
Compute time: 198.400026
Write time: 35.125583
Total time: 237.122568
09/03/2022 15:41:42     opening dtmsca was successful!
AreaDinf version 5.3.9
Input file /home/elic/fclapuyt/SedInConnect_2.3/inputs/Otemma/Otemma_filledWL2019ang.tif has projected coordinate system.
This run may take on the order of 1 minutes to complete.
This estimate is very approximate.
Run time is highly uncertain as it depends on the complexity of the input data
and speed and memory of the computer. This estimate is based on our testing on
a dual quad core Dell Xeon E5405 2.0GHz PC with 16GB RAM.
Nodata value input to create partition from file: -340282346638528859811704183484516925440.000000
Nodata value recast to float used in partition raster: -340282346638528859811704183484516925440.000000
Input file /home/elic/fclapuyt/SedInConnect_2.3/inputs/Otemma/Otemma_filledWL2019.tif has projected coordinate system.
Nodata value input to create partition from file: -99999.000000
Nodata value recast to float used in partition raster: -99999.000000
09/03/2022 15:46:42     opening acc_W was successful!
--------------------------------------------------------------------------
mpiexec was unable to find the specified executable file, and therefore
did not launch the job.  This error was first reported for process
rank 0; it may have occurred for other processes as well.

NOTE: A common cause for this error is misspelling a mpiexec command
      line parameter option (remember that mpiexec interprets the first
      unrecognized command line token as the executable).

Node:       pelican
Executable: AreaDinf
--------------------------------------------------------------------------
8 total processes failed to start
AreaDinf version 5.3.9
Input file /home/elic/fclapuyt/SedInConnect_2.3/inputs/Otemma/Otemma_filledWL2019ang.tif has projected coordinate system.
This run may take on the order of 1 minutes to complete.
This estimate is very approximate.
Run time is highly uncertain as it depends on the complexity of the input data
and speed and memory of the computer. This estimate is based on our testing on
a dual quad core Dell Xeon E5405 2.0GHz PC with 16GB RAM.
Nodata value input to create partition from file: -340282346638528859811704183484516925440.000000
Nodata value recast to float used in partition raster: -340282346638528859811704183484516925440.000000
Input file /home/elic/fclapuyt/SedInConnect_2.3/inputs/Otemma/Otemma_filledWL2019s.tif has projected coordinate system.
Nodata value input to create partition from file: -10000000000.000000
Nodata value recast to float used in partition raster: -10000000000.000000
Processors: 8
Read time: 5.898888
Compute time: 226.798826
Write time: 40.524558
Total time: 273.222272
09/03/2022 15:51:52     opening acc_S was successful!
sedinconnect_cli.py:1080: RuntimeWarning: invalid value encountered in log10
  ic_ar=numpy.log10(D_up_ar/D_down_ar)
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "sedinconnect_cli.py", line 1083, in CavalliConnectivityout
    if updown.isChecked():
NameError: name 'updown' is not defined
```

# Test Blinne

in_dtm="/home/elic/fclapuyt/SedInConnect_2.3/inputs/blinne/blinne.tif"
in_cs=2
out_w="/home/elic/fclapuyt/SedInConnect_2.3/outputs/blinne/blinne_weight.tif"
out_ic="/home/elic/fclapuyt/SedInConnect_2.3/outputs/blinne/blinne_ic.tif"

CavalliConnectivityout(in_dtm, in_cs, in_dtm, out_ic)

```bash
[Mar/09] fclapuyt@pelican 16:18:30 ~/SedInConnect_2.3 30 > python -i sedinconnect_cli.py
>>> in_dtm="/home/elic/fclapuyt/SedInConnect_2.3/inputs/blinne/blinne.tif"
>>> in_cs=2
>>> out_w="/home/elic/fclapuyt/SedInConnect_2.3/outputs/blinne/blinne_weight.tif"
>>> out_ic="/home/elic/fclapuyt/SedInConnect_2.3/outputs/blinne/blinne_ic.tif"
>>>
>>> CavalliConnectivityout(in_dtm, in_cs, in_dtm, out_ic)
09/03/2022 16:19:00     opening DTM for Connectivity was successful!
D8FlowDir version 5.3.9
Input file /home/elic/fclapuyt/SedInConnect_2.3/inputs/blinne/blinne.tif has projected coordinate system.
Nodata value input to create partition from file: -99999.000000
Nodata value recast to float used in partition raster: -99999.000000
This run may take on the order of 1 minutes to complete.
This estimate is very approximate.
Run time is highly uncertain as it depends on the complexity of the input data
and speed and memory of the computer. This estimate is based on our testing on
a dual quad core Dell Xeon E5405 2.0GHz PC with 16GB RAM.
All slopes evaluated. 0 flats to resolve.
Processors: 8
Header read time: 0.610423
Data read time: 0.150658
Compute Slope time: 0.024797
Write Slope time: 0.189779
Resolve Flat time: 0.224954
Write Flat time: 0.264076
Total time: 1.464686
09/03/2022 16:19:04     opening fd8 was successful!
DinfFlowDir version 5.3.9
Input file /home/elic/fclapuyt/SedInConnect_2.3/inputs/blinne/blinne.tif has projected coordinate system.
Nodata value input to create partition from file: -99999.000000
Nodata value recast to float used in partition raster: -99999.000000
This run may take on the order of 1 minutes to complete.
This estimate is very approximate.
Run time is highly uncertain as it depends on the complexity of the input data
and speed and memory of the computer. This estimate is based on our testing on
a dual quad core Dell Xeon E5405 2.0GHz PC with 16GB RAM.
All slopes evaluated. 0 flats to resolve.
Processors: 8
Header read time: 0.649846
Data read time: 0.148736
Compute Slope time: 0.061348
Write Slope time: 0.264761
Resolve Flat time: 0.176089
Write Flat time: 0.312638
Total time: 1.613417
09/03/2022 16:19:08     opening sd8 was successful!
09/03/2022 16:19:08     opening Weight was successful!
AreaD8 version 5.3.9
Input file /home/elic/fclapuyt/SedInConnect_2.3/inputs/blinne/blinnep.tif has projected coordinate system.
This run may take on the order of 1 minutes to complete.
This estimate is very approximate.
Run time is highly uncertain as it depends on the complexity of the input data
and speed and memory of the computer. This estimate is based on our testing on
a dual quad core Dell Xeon E5405 2.0GHz PC with 16GB RAM.
Nodata value input to create partition from file: -32768.000000
Nodata value recast to int16_t used in partition raster: -32768
Number of Processes: 8
Read time: 0.402601
Compute time: 0.402819
Write time: 0.106808
Total time: 0.912228
09/03/2022 16:19:11     Process concluded succesfully
 0.16 seconds for Weighted-Flow Length calculation with  275  iterations
AreaDinf version 5.3.9
Input file /home/elic/fclapuyt/SedInConnect_2.3/inputs/blinne/blinneang.tif has projected coordinate system.
This run may take on the order of 1 minutes to complete.
This estimate is very approximate.
Run time is highly uncertain as it depends on the complexity of the input data
and speed and memory of the computer. This estimate is based on our testing on
a dual quad core Dell Xeon E5405 2.0GHz PC with 16GB RAM.
Nodata value input to create partition from file: -340282346638528859811704183484516925440.000000
Nodata value recast to float used in partition raster: -340282346638528859811704183484516925440.000000
Processors: 8
Read time: 0.437551
Compute time: 0.364332
Write time: 0.336860
Total time: 1.138743
09/03/2022 16:19:14     opening dtmsca was successful!
AreaDinf version 5.3.9
Input file /home/elic/fclapuyt/SedInConnect_2.3/inputs/blinne/blinneang.tif has projected coordinate system.
This run may take on the order of 1 minutes to complete.
This estimate is very approximate.
Run time is highly uncertain as it depends on the complexity of the input data
and speed and memory of the computer. This estimate is based on our testing on
a dual quad core Dell Xeon E5405 2.0GHz PC with 16GB RAM.
Nodata value input to create partition from file: -340282346638528859811704183484516925440.000000
Nodata value recast to float used in partition raster: -340282346638528859811704183484516925440.000000
Input file /home/elic/fclapuyt/SedInConnect_2.3/inputs/blinne/blinne.tif has projected coordinate system.
Nodata value input to create partition from file: -99999.000000
Nodata value recast to float used in partition raster: -99999.000000
Processors: 8
Read time: 0.557960
Compute time: 0.641417
Write time: 0.309975
Total time: 1.509351
09/03/2022 16:19:18     opening acc_W was successful!
--------------------------------------------------------------------------
mpiexec was unable to find the specified executable file, and therefore
did not launch the job.  This error was first reported for process
rank 0; it may have occurred for other processes as well.

NOTE: A common cause for this error is misspelling a mpiexec command
      line parameter option (remember that mpiexec interprets the first
      unrecognized command line token as the executable).

Node:       pelican
Executable: AreaDinf
--------------------------------------------------------------------------
8 total processes failed to start
AreaDinf version 5.3.9
Input file /home/elic/fclapuyt/SedInConnect_2.3/inputs/blinne/blinneang.tif has projected coordinate system.
This run may take on the order of 1 minutes to complete.
This estimate is very approximate.
Run time is highly uncertain as it depends on the complexity of the input data
and speed and memory of the computer. This estimate is based on our testing on
a dual quad core Dell Xeon E5405 2.0GHz PC with 16GB RAM.
Nodata value input to create partition from file: -340282346638528859811704183484516925440.000000
Nodata value recast to float used in partition raster: -340282346638528859811704183484516925440.000000
Input file /home/elic/fclapuyt/SedInConnect_2.3/inputs/blinne/blinnes.tif has projected coordinate system.
Nodata value input to create partition from file: -10000000000.000000
Nodata value recast to float used in partition raster: -10000000000.000000
Processors: 8
Read time: 0.698299
Compute time: 0.472475
Write time: 0.218064
Total time: 1.388838
09/03/2022 16:19:21     opening acc_S was successful!
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "sedinconnect_cli.py", line 1083, in CavalliConnectivityout
    if updown.isChecked():
NameError: name 'updown' is not defined

```
