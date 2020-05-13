###############################################
#
#  Stutt et al. 2020
#  "A modelling framework to assess the likely effectiveness of facemasks 
# in combination with 'lock-down' in managing the COVID-19 pandemic"
#  Code for Model 1
#
###############################################


library(Rlab)
library(ggplot2)
library(RColorBrewer)

# Incubation period
# from Li et al NEJM 2020
os.par1 <- 1.434065
os.par2 <- 0.6612
onset.symptoms<-function(n) rlnorm(n, os.par1, os.par2)

# Distribution of infectiousness
#  from Xe et al Nature Medicine 2020
inf.par1<-2.115779
inf.par2<-0.6898583
inf.par3<-2.30669
infectiousness<-function(x) dgamma(x + inf.par3, inf.par1, inf.par2)

# Start generation 1 with 100 cases
num.initial.cases<-100

# Shape parameter for Negative Binomial distribution
# from Riou ans Althous 2020, Eurosurveilance
k<-0.54

# Generation when facemasks are introduced
start.wearing.face.mask<-4

# File to write results
res.file<-'results.txt'

#  Define infectiousness profile
infectiousness<-function(inc.period){
  if(inc.period>inf.par3){
  before<- pgamma(inf.par3, inf.par1, inf.par2)
  } else {
   before<- pgamma(inf.par3, inf.par1, inf.par2) -  pgamma(inf.par3 - inc.period, inf.par1, inf.par2)
  }
    after<-1-pgamma(inf.par3, inf.par1, inf.par2)
    return(c(before, after))
  }


#  Branching process model for the case when facemasks are worn ater onset of symptoms
run.model.after<-function(num.initial.cases, start.wearing.face.mask, prop.population.wearing.face.mask, effectivity.face.mask, R0, k){
  cases <- data.frame(id = 1:(num.initial.cases), exposure.time = rep(0, num.initial.cases), inc.period = onset.symptoms(num.initial.cases),   infector = 0, generation=1, new_cases = 0)
  for(generation in 2:(start.wearing.face.mask+1)){
    active.cases<-which(cases$generation==generation-1)
    for(i in 1:length(active.cases)){
      with.face.mask<-rbern(1, prop.population.wearing.face.mask)
      infectivity<- infectiousness(cases$inc.period[active.cases[i]])
      if(generation>=start.wearing.face.mask & with.face.mask==1){
         infectivity[2]<-(1-effectivity.face.mask)*infectivity[2]
      }
      secondary = rnbinom(1, size=k, mu=R0*sum(infectivity))
      if(secondary>0){
        cases$new_cases[active.cases[i]]<-secondary
        for(j in 1:secondary){
          t.new<- cases$exposure.time[active.cases[i]] + cases$inc.period[active.cases[i]] + rgamma(1, inf.par1, inf.par2) - inf.par3
          cases<-rbind(cases, data.frame(id = max(cases$id) +1, exposure.time = t.new,  inc.period = onset.symptoms(1), infector = cases$id[active.cases[i]], generation=generation, new_cases = 0))
        }
      }
    }
  }
  return(cases)
}

#  Branching process model for the case when facemasks are worn all time
run.model.before<-function(num.initial.cases, start.wearing.face.mask, prop.population.wearing.face.mask, effectivity.face.mask, R0, k){
  cases <- data.frame(id = 1:(num.initial.cases), exposure.time = rep(0, num.initial.cases), inc.period = onset.symptoms(num.initial.cases),   infector = 0, generation=1, new_cases = 0)
  for(generation in 2:(start.wearing.face.mask+1)){
    active.cases<-which(cases$generation==generation-1)
    for(i in 1:length(active.cases)){
      with.face.mask<-rbern(1, prop.population.wearing.face.mask)
      infectivity<- infectiousness(cases$inc.period[active.cases[i]])
      if(generation>=start.wearing.face.mask & with.face.mask==1){
        infectivity<-(1-effectivity.face.mask)*infectivity
      }
      secondary = rnbinom(1, size=k, mu=R0*sum(infectivity))
      if(secondary>0){
        cases$new_cases[active.cases[i]]<-secondary
        for(j in 1:secondary){
          t.new<- cases$exposure.time[active.cases[i]] + cases$inc.period[active.cases[i]] + rgamma(1, inf.par1, inf.par2) - inf.par3
          cases<-rbind(cases, data.frame(id = max(cases$id) +1, exposure.time = t.new,  inc.period = onset.symptoms(1), infector = cases$id[active.cases[i]], generation=generation, new_cases = 0))
        }
      }
    }
  }
  return(cases)
}


#  Run simulations and save output
for(i in 1:10^6){
  cat(c(i, "\t"))
  #  Sample proportion of population wearing facemask & effectiveness of facemask
  prop.population.wearing.face.mask<-sample(seq(0.05, 1, 0.1),1)
  effectivity.face.mask<-sample(seq(0.05, 1, 0.1),1)
  # Low reproduction number
  R0<-2.2 
  cases<-run.model.after(num.initial.cases, start.wearing.face.mask, prop.population.wearing.face.mask, effectivity.face.mask, R0, k)
  temp<- mean(cases$new_cases[which(cases$generation==start.wearing.face.mask)])
  cat(c(R0, 0, prop.population.wearing.face.mask, "\t", effectivity.face.mask,"\t", temp,"\n"), file = res.file, append = T, sep = "\t", eol = "\n")
  cases<-run.model.before(num.initial.cases, start.wearing.face.mask, prop.population.wearing.face.mask, effectivity.face.mask, R0, k)
  temp<- mean(cases$new_cases[which(cases$generation==start.wearing.face.mask)])
  cat(c(R0, 1, prop.population.wearing.face.mask, "\t", effectivity.face.mask,"\t", temp,"\n"), file = res.file, append = T, sep = "\t", eol = "\n")
  # High reproduction number
  R0<-4.0
  cases<-run.model.after(num.initial.cases, start.wearing.face.mask, prop.population.wearing.face.mask, effectivity.face.mask, R0, k)
  temp<- mean(cases$new_cases[which(cases$generation==start.wearing.face.mask)])
  cat(c(R0, 0, prop.population.wearing.face.mask, "\t", effectivity.face.mask,"\t", temp,"\n"), file = res.file, append = T, sep = "\t", eol = "\n")
  cases<-run.model.before(num.initial.cases, start.wearing.face.mask, prop.population.wearing.face.mask, effectivity.face.mask, R0, k)
  temp<- mean(cases$new_cases[which(cases$generation==start.wearing.face.mask)])
  cat(c(R0, 1, prop.population.wearing.face.mask, "\t", effectivity.face.mask,"\t", temp,"\n"), file = res.file, append = T, sep = "\t", eol = "\n")
}


# Upload all simulations
results<-read.table(res.file, header=F)
colnames(results)<-c("R0", "Wearing", "prop.wearing.face.mask", "effectivity.face.mask","Re")
results<-as.data.frame(results)

# Plot Figure 3
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")

grid.pnts<-unique(results[,1:4])
grid.pnts$Re<-0

for(i in 1:nrow(grid.pnts)){
  wh<-which(results[,1]==grid.pnts[i,1] & results[,2]==grid.pnts[i,2] & results[,3]==grid.pnts[i,3] & results[,4]==grid.pnts[i,4])
  grid.pnts$Re[i]<-mean(results$Re[wh])
}

grid.pnts$R0<-as.factor(grid.pnts$R0)
grid.pnts$prop.wearing.face.mask<-100*grid.pnts$prop.wearing.face.mask
grid.pnts$effectivity.face.mask<-100*grid.pnts$effectivity.face.mask


grid.pnts$Wearing[which(grid.pnts$Wearing==1)] <-"All time"
grid.pnts$Wearing[which(grid.pnts$Wearing==0)]<-"Post-symptoms"
grid.pnts$Wearing<-factor(grid.pnts$Wearing, levels=c("Post-symptoms", "All time"))


windows()
ggplot(data = grid.pnts, aes(prop.wearing.face.mask, effectivity.face.mask, fill=Re))+
  geom_tile()+
  scale_fill_gradientn(limits = c(0, max(grid.pnts$upper)), name = "Re", colours = myPalette(100)) +
  theme_minimal() +xlab("Population wearing facemask (%)") + ylab("Effectiveness of facemask (%)") +  
  facet_grid(Wearing~R0, labeller = labeller(Wearing = label_both, R0= label_both)) + theme(
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"),
    strip.text.x = element_text(size = 12),
    strip.text.y = element_text(size = 12),
    axis.text.x = element_text(size=12, angle=45),
    axis.text.y = element_text(size=12, angle=45),
    legend.title=element_text(size=14), 
    legend.text=element_text(size=12)
  ) 
