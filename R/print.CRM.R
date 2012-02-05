print.CRM <-
function(x, ...){
cat("***********************************************************************","\n")
cat("EstCRM -- An R Package for Estimating Samejima's Continuous IRT Model Parameters","\n")
cat("         ","Via Marginal Maximum Likelihood Estimation and EM Algorithm","\n")
cat("","\n")
cat("Version 1.2  2012","\n")
cat("","\n")
cat("Cengiz Zopluoglu","\n")
cat("","\n")
cat("University of Minnesota - Department of Educational Psychology","\n")
cat("","\n")
cat("zoplu001@umn.edu","\n")
cat("***********************************************************************","\n")
cat("","\n")
cat("Processing Date: ",date(),"\n")
cat("","\n")
cat("Number of Items: ",nrow(x$param),"\n")
cat("Number of Subjects: ",nrow(x$data),"\n")
cat("","\n")
cat("Item Descriptive Statistics","\n")
cat("","\n")
cat("                  ","Raw Scores","             ",
"Transformed Scores(Z scale)","\n")
cat("            ",sprintf("%6s %6s %6s %6s","Mean","SD","Min","Max"),"        ",
sprintf("%6s %6s","Mean","SD"),"\n")
for(i in 1:nrow(x$param)){
cat("     ",sprintf("%6s %6.2f %6.2f %6.2f %6.2f",colnames(x$data)[i],
x$descriptive[i,1],x$descriptive[i,2],x$descriptive[i,3],x$descriptive[i,4]),
"       ",sprintf("%6.2f %6.2f",x$descriptive[i,5],x$descriptive[i,6]),"\n")
}
cat("","\n")
cat("Iteration was terminated at EM Cycle",length(x$iterations),"\n")
cat("","\n")
cat("The Difference of Loglikelihoods Between Last Two EM Cycles=",round(x$dif,3),"\n")
cat("","\n")
cat("Largest Parameter change at the last EM Cycle is",
max(as.data.frame(x$iterations[length(x$iterations)])-as.data.frame(x$iterations[length(x$iterations)-1])),
"\n")
cat("","\n")
cat("Final Paramater Estimates:","\n")
cat("","\n")
cat("          ",sprintf("%6s %6s %8s","a","b","alpha"),"\n")
for(i in 1:nrow(x$param)){
cat("     ",sprintf("%6s %6.3f %6.3f %6.3f",colnames(x$data)[i],
x$param[i,1],x$param[i,2],x$param[i,3]),"\n")
}
}

