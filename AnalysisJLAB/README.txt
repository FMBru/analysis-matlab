Use Matlab analysis script for timing runs


Open BatchAnalyse.m in Matlab

Open batchAnalyseList.txt and enter run(s) to analyse


---- EXAMPLE for batchAnalyseList.txt ---
#run ID, Oscilloscope, DUT CH, REF CH, Detector ID (MM #), filesToAnalyse (0 is all)xxxx - add x between linesxxx  %0xxxxxxxxxxxxxxxxx
064,Pool2,4,1,2,0
x
061,Pool4,2,1,2,0
x
---- END of EXAMPLE for batchAnalyseList.txt ---

Run BatchAnalyse.m 


This will run two functions: 

First ProcessRawFiles.m 
This will loop through all files and waveforms in specified run and save extracted amp/timing information in vectors

Then AnalyseRun.m
This will produce plots and save them to Results folder





