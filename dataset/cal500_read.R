# cal500exp
x = list.files(path ="dataset/SegLabelHard")
file = NULL
for (i in 1:length(x)){
  f = read.csv(paste("dataset/SegLabelHard/", x[i], sep = ""))
  file = rbind(file, f)
}
file$Start_Time = NULL
file$End_Time = NULL
cal500exp = as.matrix(file)
save(cal500exp, file = 'cal500exp.RData')