library(readxl)
library(ggplot2)
library(ggrepel)
library(nlme)

# Scatterplots of cell number vs mtDNA area
###########

# set specific Excel coordinates of important data
base.control.row = 1
base.polg.row = 17
condition.col = 2
cellnum.col = 3

# initialise dataframe for results
new.df = data.frame()

# loop through sheets in Excel file
for(this.sheet in 2:6) {
  # sheet 2 has different coordinates for some info
  if(this.sheet == 2) {
    control.col = 4
    polg.col = 5
    base.control.col = 6
    base.polg.col = 7
  } else {
    control.col = 5
    polg.col = 6
    base.control.col = 8
    base.polg.col = 10
  }
  
  # read file and pull baseline values
  df.scat = as.data.frame(read_xlsx("Send to Iain- Scattergrams.xlsx", sheet=this.sheet))
  base.control = df.scat[base.control.row,base.control.col]
  base.polg = df.scat[base.polg.row,base.polg.col]
  base.cellnum.control = df.scat[base.control.row,cellnum.col]
  base.cellnum.polg = df.scat[base.polg.row,cellnum.col]
  
  # append baseline values to results 
  new.df = rbind(new.df, data.frame(sheet=this.sheet,condition="baseline", line="control", area=base.control, cellnum=base.cellnum.control))
  new.df = rbind(new.df, data.frame(sheet=this.sheet,condition="baseline", line="POLG", area=base.polg, cellnum=base.cellnum.polg))
  
  # identify nonzero control observations and append to results
  for(i in 1:nrow(df.scat)) {
    if(!is.na(df.scat[i,control.col])) {
      new.df = rbind(new.df, data.frame(sheet=this.sheet,condition=df.scat[i,condition.col], line="control", area=df.scat[i,control.col], cellnum=df.scat[i,cellnum.col]))
    }
  }
  
  # identify nonzero POLG observations and append to results
  for(i in 1:nrow(df.scat)) {
    if(!is.na(df.scat[i,polg.col])) {
      new.df = rbind(new.df, data.frame(sheet=this.sheet,condition=df.scat[i,condition.col], line="POLG", area=df.scat[i,polg.col], cellnum=df.scat[i,cellnum.col]))  }
  }
}

# create plot
ggplot(new.df, aes(x=cellnum, y=area, colour=line, label=condition)) + 
  geom_point(alpha=0.4) + geom_text_repel(size=2, alpha=0.4) +
  geom_point(data=new.df[new.df$condition=="baseline",], aes(x=cellnum, y=area, colour=line, label=condition), alpha=1) + 
  geom_text(data=new.df[new.df$condition=="baseline",], aes(x=cellnum, y=area, colour=line, label=condition), size=3, nudge_y = 3, alpha=1) +
  facet_wrap(~sheet) + theme_light() + theme(legend.position="none") + 
  xlab("Cell number") + ylab("mtDNA area")

# Effect of nucleoside treatment as a function of FCS concentration
########

# read Excel file
df = as.data.frame(read_xlsx("Send to Iain- ACGT conc curve Sav and HeSt.xlsx"))

# plot summary
ggplot(df, aes(x=Patcont, y=sumAreaMtDNApercell, color=factor(Concnucleosides), fill=factor(Concnucleosides))) + 
  geom_boxplot() + 
  facet_grid(FCS~Run) + theme(legend.position="none")

# reproduce Jo's summary plot
ggplot(df, aes(x=Patcont, y=sumAreaMtDNApercell, color=factor(Concnucleosides), fill=factor(Concnucleosides))) + 
  geom_bar(stat = "summary", fun = "mean", position="dodge") + 
  facet_grid(FCS~Run) + theme(legend.position="none")

# subset dataframe by FCS concentration
sub.0.1 = df[df$FCS == 0.1,]
sub.10 = df[df$FCS == 10,]
sample.0.1 = sub.0.1
sample.10 = sub.10

# linear mixed models with random effects on run number and cell line
my.lmm.0.1 = lme(log(sumAreaMtDNApercell) ~ Concnucleosides, 
                 random = list(Patcont = ~ Concnucleosides, Run = ~ Concnucleosides), 
                 data=sample.0.1, 
                 control = lmeControl(msMaxIter = 1000, msMaxEval = 1000))
my.lmm.10 = lme(log(sumAreaMtDNApercell) ~ Concnucleosides, 
                random = list(Patcont = ~ Concnucleosides, Run = ~ Concnucleosides), 
                data=sample.10, 
                control = lmeControl(msMaxIter = 1000, msMaxEval = 1000))
fit.0.1 = summary(my.lmm.0.1)
fit.10 = summary(my.lmm.10)

# pull fix effects statistics
c(fit.0.1$tTable[2,1], fit.0.1$tTable[2,2],
  fit.10$tTable[2,1], fit.10$tTable[2,2])

# effects of nucleoside combos on mtDNA area
##########

# read (big) Excel file
df = as.data.frame(read_excel("Runs470-474-476-478-484-485-488-489PGall200uM10pcFCScells_for_data_deposit.xlsx"))
# Jo's request: The key thing we are after is comparing the effect of baseline vs nucleosides for each indidivual within runs.
summary(df)

# some manual curation and subsetting
area.col = 29
colnames(df)[area.col] = "mtDNA.area"
df$Patcont[df$Patcont=="Black"] = "Control.2"
df$Patcont[df$Patcont=="Turq"] = "Control.1"
df$Patcont[df$Patcont=="HeSt"] = "POLG.1"
df$Patcont[df$Patcont=="HenSto"] = "POLG.1"
df$Patcont[df$Patcont=="SavSty"] = "POLG.2"
sub = df[c("Run", "Patcont", "Condition", "mtDNA.area")]

# summary plot
ggplot(sub, aes(x=Patcont, y=mtDNA.area, fill=Condition)) + 
  facet_wrap(~Run) + geom_boxplot()

# pull the different nucleoside combos
cond.set = unique(sub$Condition)

# for each one, try an LMM fit as above. some fail numerically; in which case catch the error
res.df = data.frame()
for(i in 2:length(cond.set)) {
  subsub = sub[sub$Condition == cond.set[1] | sub$Condition == cond.set[i],]
  result = tryCatch({ 
    this.lme = lme(log(mtDNA.area) ~ Condition, 
                   random = list(Patcont = ~ Condition, Run = ~ Condition), 
                   data=subsub, 
                   control = lmeControl(msMaxIter = 1000, msMaxEval = 1000))
    this.fit = summary(this.lme)
    res.df = rbind(res.df, data.frame(condition=cond.set[i], mean=this.fit$tTable[2,1], sd=this.fit$tTable[2,2]))
  }, error = function(err) { 
    res.df = rbind(res.df, data.frame(condition=cond.set[i], mean=NA, sd=NA))
  })
}

# compute p-values based on normal distribution of estimator
res.df$p = pnorm(0, mean=abs(res.df$mean), sd=res.df$sd)*2
