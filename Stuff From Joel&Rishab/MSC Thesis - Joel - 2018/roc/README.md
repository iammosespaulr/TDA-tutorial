# roc
Receiver Operating Characteristics<br/>
The ROC graphs are a useful tecnique for organizing classifiers and
visualizing their performance. ROC graphs are commonly used in medical
decision making.

Syntax: ROCout=roc(x,thresholds,alpha,verbose)

Input: x - This is a Nx2 data matrix. The first column is the column of the data value;
           The second column is the column of the tag: unhealthy (1) and
           healthy (0).
       Thresholds - If you want to use all unique values in x(:,1) 
           then set this variable to 0 or leave it empty; 
           else set how many unique values you want to use (min=3);
       alpha - significance level (default 0.05)
       verbose - if you want to see all reports and plots (0-no; 1-yes by
       default);

Output: if verbose = 1
        the ROCplots, the sensitivity and specificity at thresholds; the Area
        under the curve with Standard error and Confidence interval and
        comment.
        if ROCout is declared, you will have a struct:
        ROCout.AUC=Area under the curve (AUC);
        ROCout.SE=Standard error of the area;
        ROCout.ci=Confidence interval of the AUC
        ROCout.co=Cut off points
        ROCdata.xr and ROCdata.yr points for ROC plot

          Created by Giuseppe Cardillo
          giuseppe.cardillo-edta@poste.it

To cite this file, this would be an appropriate format:
Cardillo G. (2008) ROC curve: compute a Receiver Operating Characteristics curve.
http://www.mathworks.com/matlabcentral/fileexchange/19950
