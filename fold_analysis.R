#!/dors/meilerlab/apps/Linux2/x86_64/bin/Rscript

#############################################################
# Rscript for analyzing a list of folding score tables
# @author Bian Li
# Please report bugs to Bian Li via bian.li@vanderbilt.edu
#############################################################

# load the required package for parsing command line arguments
library(argparse)

# set up command line options and arguments
parser = ArgumentParser(description = "Rscript for analyzing a list of folding score tables.")
parser$add_argument("table_list", metavar = "<filename>",  nargs = "+", 
		help = "file that contains a list of folding tables to be analyzed")
parser$add_argument("--score_name", "-s", dest = "score_name", metavar = "<score_column_name>", nargs = "+",
		default = "sum", help = "name of the score based on which folding statistics is derived")
parser$add_argument("--quality_measure", "-q", dest = "quality_measure", metavar = "<quality_measure_name>", nargs = "+",
		default = "RMSD100", help = "metric used for evaluating model quality, only RMSD, RMSD100, cr6, cr12, and GDT_TS are supported")

# get command line options and arguments
args <- parser$parse_args()
table_list <- args$table_list
score <- args$score_name
measure <- args$quality_measure

# declare global variables
protein_id <- vector()
best <- vector()
top10_mean <- vector() # mean of the given measures of the top 10 models 
rmsd5 <- vector() # percentage of models with RMSD100 < 4 angstroms
rmsd8 <- vector() # percentage of models with RMSD100 < 8 angstroms
cr20 <- vector() # percentage of models with contact recovery > 20%
cr30 <- vector() # percentage of models with contact recovery > 30%
gdt_ts30 <- vector() # percentage of models with GTD_TS > 30%
gdt_ts40 <- vector() # percentage of models with GTD_TS > 40%
enrichment <- vector() # enrichment according to a specific score

# analyze folding statistics
for(files in table_list) {
	# get protein id
	con <- file(files, "r") # create file connection for reading
	tables <- readLines(con)
	for(cur_table in tables) {
		# get 4-letter PDB Code one at a time
		protein_id <- c(protein_id,substr(cur_table, 1, 4))
		# read in folding table
		score_table <- read.table(cur_table, header = TRUE)
		# remove records that have NaN
		score_table <- na.omit(score_table)
		# only interested in given quality measure
		score_table <- data.frame(score_table[, c(score, measure)])
		colnames(score_table) <- c("score", "measure")
		# sorted score table
		score_table_sorted <- data.frame()
		# variables for computing enrichment
		index <- as.integer(0.1 * nrow(score_table))
		score_sorted <- score_table[order(score_table$score, decreasing = FALSE), ]
		score_cutoff <- score_sorted$score[index]
		# collect sampling statistics
		if(measure %in% c("cr12", "cr6", "GDT_TS")) {
		  # get the best one at a time, for cr and gdt_ts, the larger the better
		  best <- c(best, max(score_table$measure, na.rm = TRUE))
		  if(measure %in% c("cr12", "cr6") ) {
		    # compute cr20 and cr30
		    cur_cr20 <- nrow(score_table[score_table$measure >= 20, ]) / nrow(score_table)
		    cr20 <- c(cr20, cur_cr20 * 100)
		    cur_cr30 <- nrow(score_table[score_table$measure >= 30, ]) / nrow(score_table)
		    cr30 <- c(cr30, cur_cr30 * 100)
		  } else if(measure == "GDT_TS") {
		    # compute gdt_ts30 and gdt_ts40
		    cur_gdt_ts30 <- nrow(score_table[score_table$measure >= 30, ]) / nrow(score_table)
		    gdt_ts30 <- c(gdt_ts30, cur_gdt_ts30 * 100)
		    cur_gdt_ts40 <- nrow(score_table[score_table$measure >= 40, ]) / nrow(score_table)
		    gdt_ts40 <- c(gdt_ts40, cur_gdt_ts40 * 100)
		  }
		  # sort measure in descending order
		  score_table_sorted <- score_table[order(score_table$measure, decreasing = TRUE), ]
		  measure_cutoff <- score_table_sorted$measure[index]
		  tp <- nrow(score_table[score_table$score <= score_cutoff & score_table$measure >= measure_cutoff, ])
		  enrichment <- c(enrichment, (tp / index) * 10)
		} else if (measure %in% c( "RMSD", "RMSD100")) {
		  # get the best one at a time, for rmsd, the smaller the better
		  best <- c(best, min(score_table$measure, na.rm = TRUE))
		  # compute rmsd5 and rmsd8
		  cur_rmsd5 <- nrow(score_table[score_table$measure <= 5, ]) / nrow(score_table)
		  rmsd5 <- c(rmsd5, cur_rmsd5 * 100)
		  cur_rmsd8 <- nrow(score_table[score_table$measure <= 8, ]) / nrow(score_table)
		  rmsd8 <- c(rmsd8, cur_rmsd8 * 100)
		  # sort rmsd in ascending order
		  score_table_sorted <- score_table[order(score_table$measure, decreasing = FALSE), ]
		  measure_cutoff <- score_table_sorted$measure[index]
		  tp <- nrow(score_table[score_table$score <= score_cutoff & score_table$measure <= measure_cutoff, ])
		  enrichment <- c(enrichment, (tp / index) * 10)
		}
		# compute the mean of top 10 models
		top10_mean <- c(top10_mean, mean(score_table_sorted$measure[1:10]))
	}
}

# print result
folding_summary <- data.frame(protein_id, best, top10_mean, enrichment)
if(measure == "cr6" || measure == "cr12") {
  folding_summary <- cbind(folding_summary, cr30, cr20)
} else if(measure == "RMSD" || measure == "RMSD100") {
  folding_summary <- cbind(folding_summary, rmsd5, rmsd8)
} else if(measure == "GDT_TS") {
  folding_summary <- cbind(folding_summary, gdt_ts40, gdt_ts30)
}
print(folding_summary)
