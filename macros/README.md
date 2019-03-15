# Description of Files

This read-me serves to provide a description of all the files in the MACROS folder. These files contain the main algorithms used in the paper, and will be useful for anyone seeking to implement the protection system.

## create_hist_lognorm.sas

Given a number of bins required, this program creates a histogram with the following outputs: bincounter; upper; lower. The histogram edges are based on a log normal distribution, with mean 11.00255 and sd 0.75275. 

This is the histogram we create for PSEO production, with numbin = 20

## create_hist_even.sas

Given a number of bins required, this program creates a histogram with the following outputs: bincounter; upper; lower. The minimum  value is 10000, and the highest value is 500,000. The top bin is 300000-500000, and the bins are evenly spaced between 10000 and 300000. 

## histogram_lognormal.sas

This program takes the following inputs: numbins, epsilon (privacy loss parameter) and an input dataset, and performs the following steps:

1. Take all the data and put them in the histogram bins, according to log-normal histogram bin constructed by create_hist_lognorm.
2. Add noise to all the counts
3. Extract the noisy percentiles of earnings, based on protected histogram counts.

## histogram_even.sas

This program takes the following inputs: numbins, epsilon (privacy loss parameter) and an input dataset, and performs the following steps:

1. Take all the data and put them in the histogram bins, according to log-normal histogram bin constructed by create_hist_even.
2. Add noise to all the counts
3. Extract the noisy percentiles of earnings, based on protected histogram counts.

## create_cdf.sas 

For the smooth sensitivity approach, we need to sample from the distribution 1/(|1+x|^4), and we have to generate that CDF analytically, which is what this program does.

## compute_ss_iteration.sas 

This program takes a list of earnings values, and calculates the smooth sensitivity for a given percentile. The required inputs are a dataset with earnings values (unsorted), the percentile to calculate SS for, and the minimum and maximum privacy loss to calculate it for. SS is then calculated between the min and max by steps of EPSSTEP, the final value to be specified.

MINEPS, MAXEPS, and EPSSTEP need to be expressed as whole numbers, since they are adjusted by dividing by 100 once in here. Therefore a value of MINEPS = 1 means that the iterations begin at epsilon=0.01.

The output is a value for the smooth sensitivity for each epsilon value.

## dp_iteration.sas 

This macro takes a dataset, variable to calculate the sensitivty on, and a list of percentiles (as well as a range of epsilon values), and outputs a dataset with protected percentile values, using the algorithm originally described in Dwork et al (2007).

