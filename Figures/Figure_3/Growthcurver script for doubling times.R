# First, load the package. 
library(growthcurver)
library(writexl)
library(ggplot2)

growthdata<- read.delim("~/UVA/Metabolic_Modeling/organized/Figures/Fig3_5X_amino_acids/Growth_Curver/growthdata5XAMINOACIDALLCFU.txt", header = TRUE)



# Load the sample growth curve data provided in the Growthcurver package.
# The first column is the time in hours, and there is one column 
# for each well in a 96-well plate.

d <- growthdata

# Let's create an output data frame to store the results in. 
# We'll create it so that it is the right size (it's faster this way!), 
# but leave it empty.
num_analyses <- length(names(d)) - 1
d_gc <- data.frame(sample = character(num_analyses),
                   k = numeric(num_analyses),
                   n0  = numeric(num_analyses),
                   r = numeric(num_analyses),
                   t_mid = numeric(num_analyses),
                   t_gen = numeric(num_analyses),
                   auc_l = numeric(num_analyses),
                   auc_e = numeric(num_analyses),
                   sigma = numeric(num_analyses),
                   stringsAsFactors = FALSE)

# Truncate or trim the input data to observations occuring in the first 20 hours.
# Remember that the times in these sample data are reported in hours. To use  
# minutes (or to trim at a different time), change the next line of code. 
# For example, if you still would like to trim at 20 hours, but your time data 
# are reported in minutes use: trim_at_time <- 20 * 60
trim_at_time <- 6*60   

# Now, loop through all of the columns in the data frame. For each column,
# run Growthcurver, save the most useful metrics in the output data frame,
# and make a plot of all the growth curve data and their best fits.

# First, create a plot for each of the wells in the 96-well plate.
# Uncomment the next line to save the plots from your 96-well plate to a 
# pdf file in the working directory.
pdf("growthcurves.pdf", height = 11, width = 3)
par(mfrow = c(10,3))
par(mar = c(0.75,0.25,0.75,0.25))
heights = c(1, 1)
widths = c(1, 1)
y_lim_max <- max(d[,setdiff(names(d), "time")])


n <- 1    # keeps track of the current row in the output data frame
for (col_name in names(d)) {
  
  # Don't process the column called "time". 
  # It contains time and not absorbance data.
  if (col_name != "time") {
    
    # Create a temporary data frame that contains just the time and current col
    d_loop <- d[, c("time", col_name)]
    
    max_value <- max(d_loop[, col_name])
   
    # Now, call Growthcurver to calculate the metrics using SummarizeGrowth
    gc_fit <- SummarizeGrowth(data_t = d_loop[, "time"], 
                              data_n = d_loop[, col_name],
                              t_trim = trim_at_time,
                              bg_correct = "none")
    
    # Now, add the metrics from this column to the next row (n) in the 
    # output data frame, and increment the row counter (n)
    d_gc$sample[n] <- col_name
    d_gc[n, 2:9] <- c(gc_fit$vals$k,
                      gc_fit$vals$n0,
                      gc_fit$vals$r,
                      gc_fit$vals$t_mid,
                      gc_fit$vals$t_gen,
                      gc_fit$vals$auc_l,
                      gc_fit$vals$auc_e,
                      gc_fit$vals$sigma)
    n <- n + 1
    
    # Finally, plot the raw data and the fitted curve

    plot(gc_fit$data$t, gc_fit$data$N, 
         pch = 20,
         ylim=c(0, max_value),
         cex = 0.6, xaxt = "n", yaxt = "n")
    text(x = trim_at_time / 4, y = max_value, labels = col_name, pos = 1)
    lines(gc_fit$data$t, predict(gc_fit$model), col = "red")
  }
}
dev.off()

write_xlsx(d_gc, "~/UVA/Metabolic_Modeling/organized/Figures/Fig3_5X_amino_acids/Growth_Curver/DoublingTimes.xlsx")


