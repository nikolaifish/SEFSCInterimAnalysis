rdat=rdat_BlackSeaBass
length_sc = length_sc
Data = new('Data')
herm = NULL
nsim=48
genus_species=NULL
Region="Southeast US Atlantic"
Fref_name = "Fmsy"
Rec="bam_recruits"
CAA_abb="all"
CAL_abb="all"
Ind_abb="all"
CV_vbK=0.001
CV_vbLinf=0.001
CV_vbt0=0.001
CV_Cat=NULL
Mat_age1_max=Mat_age1_max
length_sc=length_sc
wla_sc=NULL
wla_unit="lbs"
wla_unit_mult=1000
catch_sc=1


# Begin function

rdat <- standardize_rdat(rdat)

info <- rdat$info
parms <- rdat$parms
parm.cons <- rdat$parm.cons
parm.tvec <- rdat$parm.tvec
a.series <- rdat$a.series
t.series <- rdat$t.series
comp.mats <- rdat$comp.mats
styr <- parms$styr
endyr <- parms$endyr

Name <- gsub(" ","",str_to_title(info$species))

# MSEtool expects age-based data to begin with age 0
if(min(a.series$age)>0){
  warning(paste(Name,": Minimum age > 0. Age-based data extrapolated to age-0"))
  a.series <- data_polate(a.series,xout=0:max(a.series$age))
  a.series <- data_lim(a.series,xlim=c(0,Inf))
  a.series <- data_lim(a.series,xname=c("prop.female","prop.male","mat.female","mat.male"),xlim=c(0,1))
  a.series <- as.data.frame(a.series)
  rownames(a.series) <- a.series$age
}
age <- a.series$age

t.series <- t.series[paste(styr:endyr),]

Common_Name <- str_replace_all(Name,"(?<=[a-z])(?=[A-Z])"," ")
if(is.null(genus_species)){genus_species <- bamStockMisc[Name,"Species"]}
if(is.null(herm)){herm <- bamStockMisc[Name,"herm"]}

B.to.klb <- local({
  if(info$units.biomass%in%c("mt","metric tons")){
    out <- 2.204624}else{
      warning("units.biomass not equal to metric tons")
    }
  return(out)
})

catch.raw <- t.series$total.L.klb*catch_sc
cc.yrCat <- complete.cases(t.series$year,catch.raw) # Complete cases for catch data series
year <- t.series$year[cc.yrCat]
nyear <- length(year)
catch <- catch.raw[cc.yrCat]
recruits <- setNames(t.series$recruits,t.series$year)

# Scale BAM recruits to approximate age-0 recruits (most BAM models use age-1 for recruitment)
Nage_F0 <- expDecay(age=age,Z=a.series$M,N0=1)
bam_age_R <- min(rdat$a.series$age)
R_sc <- Nage_F0["0"]/Nage_F0[paste(bam_age_R)] # Scaling factor for recruitment
recruits_sc <- recruits*R_sc # Scaled value of BAM R0 to approximate unfished numbers at age-0

Linf <- parm.cons$Linf[8]
K <- parm.cons$K[8]
t0 <- parm.cons$t0[8]

LenCV <- parm.cons$len_cv_val[8]

M.constant <- ifelse(!is.null(parms$M.constant),
                     parms$M.constant,
                     tail(a.series$M,1))

if(Fref_name=="Fmsy"){
  Cref <- parms$msy.klb*catch_sc
  Bref <- parms$Bmsy
  Fref <- parms$Fmsy
}
if(Fref_name=="F30"){
  Cref <- parms$L.F30.klb*catch_sc
  Bref <- parms$B.F30
  Fref <- parms$F30
}

Bcurrent <- t.series[paste(endyr),"B"]
B0 <- parms$B0

SSBcurrent <- t.series[paste(endyr),"SSB"]
SSB0 <- parms$SSB0

#An estimate of absolute current vulnerable abundance, converted to klb (needs to be in same units as catch)
Abun <- sum(rdat$B.age[paste(endyr),]*B.to.klb*rdat$sel.age$sel.v.wgted.tot)*catch_sc

FMSY_M <- Fref/M.constant
BMSY_B0 <- Bref/B0
Dep <- SSBcurrent/SSB0
LHYear <- endyr

# Recruitment for years where recruitment deviations were estimated
if(!is.null(Rec)){
  if(Rec=="bam_recruits"){
    Rec <- local({
      year_nodev <- parm.tvec$year[is.na(parm.tvec$log.rec.dev)]
      recruits[paste(year_nodev)] <- NA
      matrix(data=recruits_sc,nrow=nsim,ncol=length(recruits_sc),dimnames=list("sim"=1:nsim,"year"=year))
    })
  }
}else{
  Rec <- Data@Rec
}

# Catch (Cat): Total annual catches (NOTE: DLMtool wants Cat to be a matrix)
Cat <- matrix(data=catch,nrow=nsim,ncol=length(catch),dimnames=list("sim"=1:nsim,"year"=year))

# Catch CV
if(is.null(CV_Cat)){
  CV_Cat <- matrix(0.05,nrow=nsim,ncol=nyear)
}

# Abundance Index (Ind): Relative abundance index (NOTE: DLMtool wants Ind to be a matrix)
#if(!Ind_abb[1]=="none"){
#  IndCalc <- local({
    D <- t.series[cc.yrCat,]
    x <- names(D)
    Ind_names_all <- x[grepl(pattern="U.",x=x)&grepl(pattern=".ob",x=x)]
    if(Ind_abb[1]=="all"){
      Ind_names <- Ind_names_all
    }else{
      Ind_names <- paste("U",Ind_abb,"ob",sep=".")
    }
    Ind_names <- Ind_names[Ind_names%in%names(D)] # Identify valid names
    if(length(Ind_names)==0){
      warning(paste("Ind_abb does not match any index names in the rdat t.series. Ind will be the geometric mean of all available indices:",paste(Ind_names_all,collapse=", ")))
      Ind_names <- Ind_names_all
    }
    CV_Ind_names <- paste0("cv.U.",gsub("U.|.ob","",Ind_names))
    CV_Ind_names <- CV_Ind_names[CV_Ind_names%in%names(D)] # Identify valid names

