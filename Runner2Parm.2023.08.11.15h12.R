
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~ Initial settings ~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls())
cat("\014")
set.seed(1234)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(stargazer)
source("Estimation.R")
source("GoF_Pareto_MME.R")
source("DistsForSim.2023.08.11.15h42.R")

starttime <- format(Sys.time(), "%Y%m%d-%Hh%Mm%Ss")

powersWS <- function(TS,alpha){
  TS_X = TS[,1]
  TS_S = TS[,2]
  cv   = quantile(TS_S,1-alpha)
  out  = mean(TS_X>cv)
  return(out)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

n       = 20
Range   = c(73:104) #c(3,9:12,13:16,25:28,29:32,33:36,46:54,64:72)  #distributions used
MC      = 2e4
alpha   = 0.05
meth    = "MLE" #or "MME"
m.vals  = sort(unique(c(3,4))) #sort(unique(c(2,3,4,(1:9)*n/10)))
num.m   = length(m.vals)
a.vals  = 0.5#c(0.25,0.5,0.75,1,1.5,2,5)
num.a   = length(a.vals) 
ntests0 = 13
ntests  = ntests0 + 2*num.m*num.a
Results = matrix(0,length(Range),ntests + 1)  # (0, number of distributions, number of tests + 1)
nm      = numeric(length(Range))
tme1    = proc.time()
#path   <- getwd()
#filnam <- paste(path,"/InverseGaussian_n=",n,"_MC=",MC,"_Dist=",Range[1],"to",Range[length(Range)],".txt",sep="")
#cat("\n",file=filnam,append=TRUE)

count = 0
pb    = txtProgressBar(min=0,max=length(Range),style=3)
DistNumcntr <- 1
for (DistNum in Range){#START FOR LOOP (DistNum in Range)
  p       = numeric(ntests) 
  T1      = array(0,c(MC,2,num.m*num.a))
  T2      = array(0,c(MC,2,num.m*num.a))
  TS_KS   = 
    TS_CV   = 
    # TS_AD   = 
    TS_MA   = 
    TS_ZA   = 
    TS_ZB   = 
    # TS_ZC   = 
    TS_KL_1 = 
    TS_KL_2 = 
    TS_DK   = 
    TS_ME_1 = 
    TS_ME_2 = 
    TS_M2_1 = 
    TS_M2_2 = 
    TS_A1   = matrix(0,MC,2)
  
  
  for (j in 1:MC){#START FOR LOOP (j in 1:MC)
    simdata     = SimFromDist(DistNum,n)
    X           = sort(simdata$X)
    nm[DistNumcntr] = simdata$nm
    # betaH       = ifelse(meth=="MME",ParetoMME1(X),ParetoMLE1(X)) # MME or MLE of beta is calculated
    betaH       = ifelse(meth=="MME",ParetoMME2(X)$betaH  , ParetoMLE2(X)$betaH) # MME or MLE of beta is calculated
    sigmaH      = ifelse(meth=="MME",ParetoMME2(X)$sigmaH ,ParetoMLE2(X)$sigmaH) # MME or MLE of sigma is calculated
    Xstar       = sort(rpareto(n,betaH))
    # betaHS      = ifelse(meth=="MME",ParetoMME1(Xstar),ParetoMLE1(Xstar)) # MME or MLE of beta is calculated
    betaHS      = ifelse(meth=="MME",ParetoMME2(Xstar)$betaH ,ParetoMLE2(Xstar)$betaH ) # MME or MLE of beta is calculated
    sigmaHS     = ifelse(meth=="MME",ParetoMME2(Xstar)$sigmaH,ParetoMLE2(Xstar)$sigmaH) # MME or MLE of sigma is calculated
    TS_KS[j,]   = c(KSn(X/sigmaH,betaH),    KSn(Xstar/sigmaHS,betaHS))
    TS_CV[j,]   = c(CVn(X/sigmaH,betaH),    CVn(Xstar/sigmaHS,betaHS))
    # TS_AD[j,]   = c(ADn(X/sigmaH,betaH),    ADn(Xstar/sigmaHS,betaHS))
    TS_MA[j,]   = c(MAn(X/sigmaH,betaH),    MAn(Xstar/sigmaHS,betaHS))
    TS_ZA[j,]   = c(ZAn(X/sigmaH,betaH),    ZAn(Xstar/sigmaHS,betaHS))
    TS_ZB[j,]   = c(ZBn(X/sigmaH,betaH),    ZBn(Xstar/sigmaHS,betaHS))
    # TS_ZC[j,]   = c(ZCn(X/sigmaH,betaH),    ZCn(Xstar/sigmaHS,betaHS))
    TS_KL_1[j,] = c(KLn(X/sigmaH,1,betaH),  KLn(Xstar/sigmaHS,1,betaHS))
    TS_KL_2[j,] = c(KLn(X/sigmaH,10,betaH), KLn(Xstar/sigmaHS,10,betaHS))
    TS_DK[j,]   = c(DKn(X/sigmaH,betaH),    DKn(Xstar/sigmaHS,betaHS))
    TS_ME_1[j,] = c(MEn(X/sigmaH,0.5,betaH),MEn(Xstar/sigmaHS,0.5,betaHS))
    TS_ME_2[j,] = c(MEn(X/sigmaH,1,betaH),  MEn(Xstar/sigmaHS,1,betaHS))
    TS_M2_1[j,] = c(M2n(X/sigmaH,0.5,betaH),M2n(Xstar/sigmaHS,0.5,betaHS))
    TS_M2_2[j,] = c(M2n(X/sigmaH,2,betaH),  M2n(Xstar/sigmaHS,2,betaHS))
    TS_A1[j,]   = c(A1n(X/sigmaH,2),        A1n(Xstar/sigmaHS,2))
    tmpcntr <- 1
    for(mm in m.vals){ #START FOR LOOP (mm in m.vals)
      for(aa in a.vals){ #START FOR LOOP (aa in a.vals)
        T1[j,,tmpcntr] = c(Tnm1(X/sigmaH,mm,aa,betaH),Tnm1(Xstar/sigmaHS,mm,aa,betaHS))
        T2[j,,tmpcntr] = c(Tnm2(X/sigmaH,mm,aa,betaH),Tnm2(Xstar/sigmaHS,mm,aa,betaHS))
        tmpcntr <- tmpcntr + 1
      }#END FOR LOOP (aa in a.vals)
    }#END FOR LOOP (mm in m.vals)
  }#END FOR LOOP (j in 1:MC)
  
  # p[1]  = powersWS(TS_KS,alpha)
  # p[2]  = powersWS(TS_CV,alpha)
  # p[3]  = powersWS(TS_AD,alpha)
  # p[4]  = powersWS(TS_MA,alpha)
  # p[5]  = powersWS(TS_ZA,alpha)
  # p[6]  = powersWS(TS_ZB,alpha)
  # p[7]  = powersWS(TS_ZC,alpha)
  # p[8]  = powersWS(TS_KL_1,alpha)
  # p[9]  = powersWS(TS_KL_2,alpha)
  # p[10] = powersWS(TS_DK,alpha)
  # p[11] = powersWS(TS_ME_1,alpha)
  # p[12] = powersWS(TS_ME_2,alpha)
  # p[13] = powersWS(TS_M2_1,alpha)
  # p[14] = powersWS(TS_M2_2,alpha)
  # p[15] = powersWS(TS_A1,alpha)
  
  p[1]  = powersWS(TS_KS,alpha)
  p[2]  = powersWS(TS_CV,alpha)
  p[3]  = powersWS(TS_MA,alpha)
  p[4]  = powersWS(TS_ZA,alpha)
  p[5]  = powersWS(TS_ZB,alpha)
  p[6]  = powersWS(TS_KL_1,alpha)
  p[7]  = powersWS(TS_KL_2,alpha)
  p[8] = powersWS(TS_DK,alpha)
  p[9] = powersWS(TS_ME_1,alpha)
  p[10] = powersWS(TS_ME_2,alpha)
  p[11] = powersWS(TS_M2_1,alpha)
  p[12] = powersWS(TS_M2_2,alpha)
  p[13] = powersWS(TS_A1,alpha)
  tmpcntr <- ntests0 + 1 #13
  for(mm in m.vals){ #START FOR LOOP (mm in m.vals)
    for(aa in a.vals){ #START FOR LOOP (aa in a.vals)
      p[tmpcntr]                 = powersWS(T1[,,tmpcntr - ntests0],alpha)
      p[tmpcntr+(num.m*num.a)]   = powersWS(T2[,,tmpcntr - ntests0],alpha)
      tmpcntr                    = tmpcntr + 1
    }#END FOR LOOP (aa in a.vals)
  }#END FOR LOOP (mm in m.vals)
  
  Results[DistNumcntr,] = c(DistNum,p)
  #	cat(Results[DistNum,],"\n",file=filnam,append=TRUE)
  count = count+1
  imagename <- paste("Pareto_n=",n,"_MC=",MC,"_parm=1_meth=",meth,".RData",sep="")
  save.image(imagename)
  DistNumcntr <- DistNumcntr + 1
  setTxtProgressBar(pb,count)
}#END FOR LOOP (DistNum in Range)

tme2 = proc.time()
tme  = tme2[3]-tme1[3]


# statnm <- c("Distribution",
#             "TS_KS",
#             "TS_CV",
#             "TS_AD",
#             "TS_MA",
#             "TS_ZA",
#             "TS_ZB",
#             "TS_ZC",
#             "TS_KL_1",
#             "TS_KL_2",
#             "TS_DK",
#             "TS_ME_1",
#             "TS_ME_2",
#             "TS_M2_1",
#             "TS_M2_2",
#             "TS_A1",
#             paste0("T1_",paste0("m=",apply(expand.grid(m.vals,a.vals),1,paste0,collapse=",a="))),
#             paste0("T2_",paste0("m=",apply(expand.grid(m.vals,a.vals),1,paste0,collapse=",a="))))


statnm <- c("Distribution",
            "TS_KS",
            "TS_CV",
            "TS_MA",
            "TS_ZA",
            "TS_ZB",
            "TS_KL_1",
            "TS_KL_2",
            "TS_DK",
            "TS_ME_1",
            "TS_ME_2",
            "TS_M2_1",
            "TS_M2_2",
            "TS_A1",
            paste0("T1_",paste0("m=",apply(expand.grid(m.vals,a.vals),1,paste0,collapse=",a="))),
            paste0("T2_",paste0("m=",apply(expand.grid(m.vals,a.vals),1,paste0,collapse=",a="))))


out     = round(Results*100)
out[,1] = out[,1]/100
dimnames(out) <- list(nm,statnm)
# out

# stargazer(out)
# out

imagename <- paste("Pareto_n=",n,"_MC=",MC,"_parm=1_meth=",meth,".RData",sep="")
save.image(imagename)
# load("ParetoExisting_n=20_MC=10.RData")

#tme/MC*50*1000/60/60

# out[Range,]

write.csv(out,file=paste0(starttime,"_Pareto_n=",n,"_MC=",MC,"_parm=1_meth=",meth,"_HA=",min(Range),"-",max(Range),".csv"))

cat("\nTime elapsed:\t",round(tme/60,2),"minutes")



