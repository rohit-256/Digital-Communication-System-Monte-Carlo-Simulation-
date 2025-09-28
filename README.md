main1.m simulates MQAM and MPSK systems over AWGN channels.
main1(constellation,M) plots the analytical and empirical SER vs SNR plot for a communication system using an M point constellation.
constellation = 'MPSK' or 'MQAM'

######################################################################################################################################
main2 simulates a communication system with a custom 2_D constellation specified by a csv file.
eg: main2(<filename.csv>)

testData.csv is given as a sample for the format
