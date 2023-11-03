library(readxl)
library(ggplot2)
library(nlme)

###########

base.control.row = 1
base.polg.row = 17
condition.col = 2

new.df = data.frame()

for(this.sheet in 2:5) {
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
df.scat = as.data.frame(read_xlsx("Send to Iain- Scattergrams.xlsx", sheet=this.sheet))
base.control = df.scat[base.control.row,base.control.col]
base.polg = df.scat[base.polg.row,base.polg.col]
for(i in 1:nrow(df.scat)) {
  if(!is.na(df.scat[i,control.col])) {
    new.df = rbind(new.df, data.frame(sheet=this.sheet,condition=df.scat[i,condition.col], control=df.scat[i,control.col]/base.control, polg=0))
  }
}
for(i in 1:nrow(df.scat)) {
  if(!is.na(df.scat[i,polg.col])) {
    ref = which(new.df$condition == df.scat[i,condition.col] & new.df$sheet==this.sheet)
    new.df$polg[ref] = df.scat[i,polg.col]/base.polg
  }
}
}
ggplot(new.df, aes(x=control, y=polg, label=condition)) + geom_point() + geom_text() + facet_wrap(~sheet)

########

df = as.data.frame(read_xlsx("Send to Iain- ACGT conc curve Sav and HeSt.xlsx"))

ggplot(df, aes(x=Patcont, y=sumAreaMtDNApercell, color=factor(Concnucleosides), fill=factor(Concnucleosides))) + 
  geom_boxplot() + 
  facet_grid(FCS~Run) + theme(legend.position="none")

# reproduce Jo's summary plot
ggplot(df, aes(x=Patcont, y=sumAreaMtDNApercell, color=factor(Concnucleosides), fill=factor(Concnucleosides))) + 
geom_bar(stat = "summary", fun = "mean", position="dodge") + 
  facet_grid(FCS~Run) + theme(legend.position="none")

set.seed(1)
sub.0.1 = df[df$FCS == 0.1,]
sub.10 = df[df$FCS == 10,]
sample.0.1 = sub.0.1[sample(1:nrow(sub.0.1), 1000),]
sample.10 = sub.10[sample(1:nrow(sub.10), 1000),]
sample.0.1 = sub.0.1
sample.10 = sub.10

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

c(fit.0.1$tTable[2,1], fit.0.1$tTable[2,2],
  fit.10$tTable[2,1], fit.10$tTable[2,2])

##########

df = as.data.frame(read_excel("Runs470-474-476-478-484-485-488-489PGall200uM10pcFCScells_for_data_deposit.xlsx"))
#The key thing we are after is comparing the effect of baseline vs nucleosides for each indidivual within runs.
summary(df)
area.col = 29
colnames(df)[area.col] = "mtDNA.area"
df$Patcont[df$Patcont=="Black"] = "Control.2"
df$Patcont[df$Patcont=="Turq"] = "Control.1"
df$Patcont[df$Patcont=="HeSt"] = "POLG.1"
df$Patcont[df$Patcont=="HenSto"] = "POLG.1"
df$Patcont[df$Patcont=="SavSty"] = "POLG.2"
sub = df[c("Run", "Patcont", "Condition", "mtDNA.area")]

ggplot(sub, aes(x=Patcont, y=mtDNA.area, fill=Condition)) + 
  facet_wrap(~Run) + geom_boxplot()

cond.set = unique(sub$Condition)

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


