################################################################################
# This file replicates Tables 4 and 5 from Diebold and Shin (2019, IJF):
## First version: 12/1/2019
## Current version: 3/25/2024
################################################################################
# 1. Directly motivated by the ex-post peLASSO solution, we select a small number (N) of the "best" forecasts and average them. Inspired by this, we propose the 'Best Nmax Average' approach. This approach searches over all possible combinations of K forecasters where N is less than or equal to Nmax. The appeal of this approach lies in the fact that, unlike the two-step peLASSO approach, we do not have to choose a tuning parameter. To replicate Table 4, set `table_to_replicate = 4`.

# 2. We investigate the performance of the same method that adaptively selects the estimation window size (as shown in Table 5). In our application, this approach enhances the forecasting performance. However, as this approach searches over the set of window sizes, it requires a longer computation time.

# 3. The common part of the table (RMSEs from simple average, best, worst etc forecasters) can be computed from the following file: "replication_DS2019_tab_common.R"

# If you have any questions, please contact us:
# Francis X. Diebold (fdiebold@upenn.edu)
# Minchul Shin (visiblehand@gmail.com)
################################################################################

# Housekeeping
library(Matrix);
library(pracma);
rm(list = ls());

# Main option for the file
table_to_replicate = 4; #4 for table 4; 5 for table 5;

if (table_to_replicate == 4){
  wgrid = 20; #our default in the paper
} else if (table_to_replicate == 5){
  wgrid = cbind(1,2,3,4,5,6,7,8,9,10,15,20,25,30,35,40); #grid for window size
} else {
  stop("Only table_to_replicate = {4, 5} is supported.")
}


######################################################################
## Setting
tevalini = 6; # starting point of prediction evaluation
tevalini0 = 2; # starting point of prediction evaluation, before the window selection

kmax = 6; # the maximum number of forecasters included in the pool

nw = length(wgrid); # number of windows under consideration

## read data
loadfilename = "H1_gdp.csv";
savefilename = "bestavg_W_H1_gdp.Rdata";

GDP = read.csv(file = loadfilename, header = TRUE);
GDP = subset(GDP, select = -c(Date.of.forecast) );
GDP = as.matrix(GDP);
row = nrow(GDP);
column = ncol(GDP);

response0 = as.matrix(GDP[, column]); #last column is the actual data
input0    = as.matrix(GDP[, 2:(column-1)]); # rest are the forecasts

######################################################################
## Actual program starts here
bestcombmaxW.n.t = matrix(NA, nrow=row, ncol=nw);
bestcombmaxW.avg.t = matrix(NA, nrow=row, ncol=nw);
bestcombmaxW.err.t = matrix(NA, nrow=row, ncol=nw);

## Loop over window size
for (wind in 1:nw){

  w = wgrid[wind];
  print(cbind("Window: ", wind, " out of ", nw));

  ###########################
  ## Part 2-1) best performing groups (N choose K)
  bestcomb.id.t = vector("list", kmax); # initialize list
  for (k in 1:kmax){
    bestcomb.id.t[[k]] = matrix(data= NA, nrow = row, ncol= k);
  }

  bestcomb.err.t = matrix(data = NA, nrow = row, ncol = kmax);
  bestcomb.rmse  = matrix(NA,nrow=kmax,ncol=1);
  bestset.rmse.t = matrix(data = NA, nrow = row, ncol = 6); #historical best

  ## Actual iteration over k
  for (k in 1:kmax){

    print(cbind(k," is processing out of kmax = ", kmax));

    # actual iteration over t
    t = 0;
    for (t in tevalini0:row){

      # set-up
      t0 = max(1, t-w);
      input = input0[t0:(t-1), , drop=F];
      response = response0[t0:(t-1), , drop=F];

      # find the best combination from the past
      temp_id   = combn(nid,k);
      temp_n    = dim(temp_id)[2];
      temp_rmse = matrix(data = NA, nrow = temp_n, ncol=1);
      for (i  in 1:temp_n){
        temp_rmse[i] = sqrt(mean( (response-rowMeans(input[, temp_id[,i], drop=F]))^2));
      }

      # compute average based on the best combination
      bestset.rmse.t[t,k] = min(temp_rmse); # best based on the historical data (t-w, t-1)
      temp_best_id = temp_id[, which.min(temp_rmse)];
      temp_best_avg = mean( input0[t, temp_best_id]);
      temp_best_err = response0[t, ] - temp_best_avg; # forecast error at t (of best combination)

      # store
      bestcomb.err.t[t,k] = temp_best_err;
      bestcomb.id.t[[k]][t, ] = temp_best_id;

    } # loop over t
    bestcomb.rmse[k] = sqrt(mean(bestcomb.err.t[,k] ^2, na.rm = TRUE));

  } # loop over k

  ###########################
  ## Part 2, <Nmax approach (only computes kmax case)
  # search for best k <= kmax based on the historical performance at time t (to record error at t)
  t = 0;
  for (t in tevalini0:row){

    # set-up
    t0 = max(1, t-w);
    input = input0[t0:(t-1), , drop=F];
    response = response0[t0:(t-1), , drop=F];

    # historical best combinations (based on data up to t-1)
    # - we've already recorded best comb for each k; we only need to compare them
    temp_rmse2 = matrix(NA,nrow=kmax,ncol=1);
    for (k in 1:kmax){
      temp_rmse2[k] = sqrt( mean( (response - rowMeans(input[,bestcomb.id.t[[k]][t, ,drop=F],drop=F],na.rm=T))^2 ) );
    }

    # best Nmax id at t with w
    temp_id2 = which.min(temp_rmse2[1:kmax]);
    temp_best_avg = mean(input0[t,bestcomb.id.t[[temp_id2]][t,]]);
    temp_best_err = response0[t,] - temp_best_avg;

    # store
    bestcombmaxW.n.t[t,wind]   = temp_id2;
    bestcombmaxW.avg.t[t,wind] = temp_best_avg;
    bestcombmaxW.err.t[t,wind] = temp_best_err; # error from bestavg_kmax approach at t (that is, build everything using data up to t-1; score it using realization at t)
    } # loop over t
} # loop over w
sqrt(colMeans(bestcombmaxW.err.t[tevalini:row, ,drop=F] ^ 2, na.rm = TRUE));

#########################################################
## Part 3: Window size selection
qgrid = 90; # we set this large enough that we include all data observations for the window size selection
#qgrid = wgrid; # alternatively
nq = length(qgrid);

bestcombmaxWW.err.t = matrix(NA,nrow =row,ncol=nq);
bestcombmaxWW.n.t   = matrix(NA,nrow =row,ncol=nq);
bestcombmaxWW.w.t   = matrix(NA,nrow =row,ncol=nq);

for (qind in 1:nq){
  q = qgrid[qind];

  tevaliniW = 6;
  for (t in tevaliniW:row){

    # window selection based on historical performance
    t0 = max(1,t-q); #recursive
    t1 = t-1;
    temp_h_err  = bestcombmaxW.err.t[t0:t1, , drop=F]; #historical error
    temp_h_rmse = sqrt(colMeans(temp_h_err^2 ,na.rm=T)) #historical rmse

    # prediction for tomorrow
    temp_best_w = which.min(temp_h_rmse);
    bestcombmaxWW.err.t[t,qind] = bestcombmaxW.err.t[t,temp_best_w];
    bestcombmaxWW.n.t[t,qind]   = bestcombmaxW.n.t[t,temp_best_w];
    bestcombmaxWW.w.t[t,qind]   = wgrid[temp_best_w];
  }
}
# performance of adaptively chosen W
bestcombmaxWW.rmse = sqrt(colMeans(bestcombmaxWW.err.t[tevalini:row, ,drop=F] ^ 2, na.rm = TRUE));
# performance of all w
bestcombmaxW.rmse = sqrt(colMeans(bestcombmaxW.err.t[tevalini:row, ,drop=F] ^ 2, na.rm = TRUE));

## Save
save.image(file = savefilename);

# ==============================================================================
print("final results");
print(bestcombmaxW.rmse);


















