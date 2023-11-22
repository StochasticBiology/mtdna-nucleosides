library(readxl)
library(ggplot2)

# read data
df = as.data.frame(read_excel("Fratter copy no supplements Nov 2020 for SPSS.xlsx"))

# day 0 "ddC" serves as the baseline for all measurements
# "DMSO" is control and "ATGC" is treatment
# "Plate" gives the index of independent replicates

# set baseline references
base.ATGC = "ddC"
base.DMSO = "ddC"

# initialise normalised mtDNA variable
df$normed = NA
# normalise all appropriate rows in the dataframe by baseline reference
for(i in 1:nrow(df)) {
  this.ref = 0
  if(df$cond[i]=="ATGC" | (df$cond[i]==base.ATGC & df$day[i]==0)) {
    this.ref = mean(df$mtDNA[  df$patcont==df$patcont[i] &
                               df$cond==base.ATGC &
                               df$day == 0])
  } else if(df$cond[i]=="DMSO" | (df$cond[i]==base.DMSO & df$day[i]==0)) {
    this.ref = mean(df$mtDNA[  df$patcont==df$patcont[i] &
                               df$cond==base.DMSO &
                               df$day == 0])
  }
  if(this.ref == 0) {
    df$normed[i] = NA
  } else {
    df$normed[i] = df$mtDNA[i]/this.ref
  }
}

# visualise raw data
ggplot(df, aes(x=day, y=mtDNA, color=cond)) + 
  geom_point() + 
  facet_wrap(~patcont)

# visualise and save normalised data
sf = 2
png("time-series.png", width=500*sf, height=300*sf, res=72*sf)
ggplot(df[df$cond != "NTC",], aes(x=day, y=normed, color=cond)) + 
  geom_point(size=1.5) + xlab("Day") + ylab("mtDNA normalised by Day 0 mean") +
  facet_wrap(~patcont) + theme_light()
dev.off()

# ANOVAs blocking "day" and reporting effect of treatment on normalised mtDNA
# (for each cell line)
# not the best -- time series measurements aren't independent
res = list()
for(patcont in unique(df$patcont)) {
  aov.control = aov(normed ~ day + cond, 
    data=df[df$patcont==patcont & df$cond %in% c("DMSO", "ATGC"),])
  res[[length(res)+1]] = summary(aov.control)
}
unique(df$patcont)

# individual t-tests for ATGC vs DMSO at each timepoint
# (for each cell line)
# low power (n=2 for each)
res.df = data.frame()
for(patcont in unique(df$patcont)) {
 sub = df[df$patcont==patcont & df$cond != "NTC",]
 for(day in c(14,21)) {
   my.t = t.test(sub$normed[sub$day == day & sub$cond=="ATGC"],
          sub$normed[sub$day == day & sub$cond=="DMSO"])
   res.df = rbind(res.df, data.frame(patcont=patcont,day=day, p=my.t$p.value))
 }
}

