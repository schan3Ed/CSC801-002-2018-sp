


library(tableHTML)
library(xtable, help, pos = 2, lib.loc = NULL)

df <- read.table("/home/chrx/Documents/git/CSC801-002-2018-sp/_course_lib_r/data/fg-asym-pi_needles_test1.txt")
#print (df)
tableHTML(df)
write_tableHTML(tableHTML(df), file = "test.html")

#print(xtable(df), type = "html")