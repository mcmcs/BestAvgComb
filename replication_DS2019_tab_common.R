################################################################################
# This file replicates the common part of the tables from Diebold and Shin (2019, IJF):
# See replication_DS2019_tab045.R for replicating the main part of the paper.
# First version: 12/1/2019
# Current version: 3/25/2024
################################################################################
# If you have any questions, please contact us:
# Francis X. Diebold (fdiebold@upenn.edu)
# Minchul Shin (visiblehand@gmail.com)
################################################################################

## House keeping
rm(list = ls())

## Packages
library(Matrix);
library(pracma);

################################################################################
## Setting parameters
w = 20; # window size
tevalini  = 6; # starting point of prediction evaluation (actual evaluation)
tevalini0 = 5; # starting point of prediction exercise (actual prediction)
tevalini1 = 6; #initial evaluation period

ugrid = exp(seq(15, -15, length = 200));
ugrid2 = exp(seq(15, -15, length = 200)); #two-step method grid
## this is what we use in paper

################################################################################
## read data
loadfilename = "H1_gdp.csv";

GDP = read.csv(file = loadfilename, header = TRUE);
GDP = subset(GDP, select = -c(Date.of.forecast) );

GDP = as.matrix(GDP);
row = nrow(GDP);
column = ncol(GDP);

input0 = as.matrix(GDP[, 2:(column - 1)]);
response0 = as.matrix(GDP[, column]);

nid = ncol(input0);
id = seq(1, nid);

input = matrix(data = NA);
response = as.vector(NULL);

################################################################################
##Equal Weights
################################################################################

## Preparation (common part of the table)
# simple average
eq.err.t = matrix(data = NA, nrow = row, ncol = 1);
for (t in tevalini1:row) {
  eq.err.t[t,] = GDP[t, column]-mean(GDP[t, 2:(column - 1)]);
};
eq.rmse = sqrt(mean(eq.err.t ^ 2, na.rm=T));

# individual forecasters
ind.err.t = matrix(data = NA, nrow = row, ncol = nid);
for (t in tevalini1:row) {
  for (i in 1:nid) {
    ind.err.t[t, i] = response0[t, 1] - input0[t, i];
  };
};
ind.rmse = sqrt(colMeans(ind.err.t ^ 2, na.rm = T));

# Create a table
tab_ind_c0 = rbind("Best Individual", "90% Individual", "Median (50%) Individual", "10% Individual", "Worst Individual");
tab_ind_c1 = matrix(quantile(ind.rmse, cbind(0,0.1,0.5,0.9,1)), nrow = 5);
tab_ind_c2 = matrix(cbind(1,1,1,1,1), nrow=5, ncol=1);

tab_eq_c0 = "Simple Average";
tab_eq_c1 = eq.rmse;
tab_eq_c2 = nid;

tab_head   = cbind("Comparison Group", "RMSE", "#");
tab_common = as.data.frame.matrix(rbind(tab_head, cbind(tab_ind_c0, tab_ind_c1, tab_ind_c2), cbind(tab_eq_c0, tab_eq_c1, tab_eq_c2)));

# print table to the screen
print(tab_common)


