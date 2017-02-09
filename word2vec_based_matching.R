# load data and pre-processing
load("CCR_APDC.Rdata")
library(dplyr)
library(lubridate)
library(ggplot2)
# load information and ICD codes
dt <- CCR_APDC %>% filter(acuteflg=="Yes" & !substr(csepmode,1,4)=="Died") %>%
  select(PPN,LHD_hosp,YearAdmit,agegrp,readmit,sex,c(86:164)) %>% 
  droplevels()
# load dictionary
dict <- read.csv("dict.csv", header=T, sep=",")
# convert to ICD codes
icd_col <- 7
for(i in icd_col:ncol(dt))
{
  dt[,i] <- dict$code_id[match(dt[,i], dict$ascii_desc)]
}
saveRDS(dt, "original_dt.rds")
# count frequency
library(dplyr)
cnt <- dt[,c(icd_col:ncol(dt))] %>% as.matrix %>% table() %>% as.data.frame() %>% arrange(desc(Freq))
saveRDS(cnt, "freq.rds")
# remove ICD codes less than 30
cnt <- cnt[which(cnt$Freq>=30),]
dt.icd <- apply(dt[,c(icd_col:ncol(dt))], 2, function(x) {
  cnt[,1][match(x,cnt[,1])]
})
dt.icd <- as.data.frame(dt.icd)
dt <- cbind(dt[,c(1:(icd_col-1))],dt.icd)
# remove NA by shifting cells
# x is a combination of not NA and NA
# we need to transpose
dt <- t(apply(dt, MARGIN=1, FUN=function(x) {
  x <- c(x[!is.na(x)],x[is.na(x)])
}))
dt <- as.data.frame(dt)
# remove empty columns
dt <- dt[,colSums(!is.na(dt))>0]
# get # icd codes in each patient
size <- rowSums(!is.na(dt))
hist(size, breaks=100)
# save data
saveRDS(dt, "freq_30_dt.rds")
# save new data
dt_new <- CCR_APDC %>% filter(acuteflg=="Yes" & !substr(csepmode,1,4)=="Died") %>%
  select(sepdate,dth_date) %>% droplevels()
saveRDS(dt_new, "new_dt.rds")
# load data
dt <- readRDS("freq_30_dt.rds")
dt_new <- readRDS("new_dt.rds")
dt <- cbind(dt_new,dt)
# get visits having at least 1 icd codes
icd_col <- 9
df <- dt[rowSums(!is.na(dt))>=icd_col,]
df[,1] <- as.Date(df[,1], format = "%d/%m/%Y")
df[,2] <- as.Date(df[,2], format = "%d/%m/%Y")
cs_date <- as.Date("31/12/2008", format = "%d/%m/%Y")
# cut-off date of cancer study
df <- df %>% filter(sepdate<=cs_date)
# rename variables
df <- df %>% rename(pid=V1,hos=V2,year=V3,age=V4,readm=V5,gender=V6)
# create file for word2vec
write.table(df[,c(icd_col:ncol(df))], "input.txt", quote=F, sep=" ", row.names=F, col.names=F, na="")

# analyze original data
dt <- readRDS("original_dt.rds")
# count admissions
nrow(dt)
# count unique patients
length(unique(dt$PPN))
# statistic of year admit
summary(dt$YearAdmit)
cnt <- readRDS("freq.rds")
# count unique icd codes
nrow(cnt)
# statistic of icd codes
summary(cnt$Freq)
# analyze subset data for control matching
# count admissions
nrow(df)
# count unique patients
length(unique(df$pid))
# statistic of ICD codes in admissions
icd_size <- rowSums(!is.na(df[,c(icd_col:ncol(df))]))
hist(size, breaks=100)
summary(size)
blank_theme <- theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid=element_blank(),
    axis.ticks.length = unit(.5, "cm"),
    axis.text = element_text(size = 40)
  )
ggplot() + aes(icd_size) + geom_histogram(binwidth = 0.5) + blank_theme
  #xlab("# of ICD codes") + ylab("# of admissions")
# readmission status
table(df$readm)*100/nrow(df)
readm_dt <- data.frame(
  readmission_status = c("Missing (0.01%)","Readmitted within 28 days to another facility (4.07%)",
            "Readmitted within 28 days to the same facility (23.03%)",
            "Not formally readmitted within 28 days (72.89%)"),
  value = c(0.01,4.07,23.03,72.89)
)
readm_dt %>% ggplot(aes(x="", y=value, fill=readmission_status)) +
  geom_bar(width = 1, stat = "identity") + coord_polar("y", start=0) +
  scale_fill_brewer(palette="Dark2") + theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_text(size = 40),
    #legend.justification = c(1, 1), 
    legend.position = "none"
  )
# length of survival
df_death <- df %>% filter(!is.na(dth_date))
df_death$LOS <- as.numeric(as.duration(df_death[,"dth_date"]-df_death[,"sepdate"])/(3600*24*365))
ggplot() + aes(df_death$LOS[df_death$LOS>0]) + geom_histogram(binwidth = 0.05) + blank_theme
  #xlab("length of survival (years)") + ylab("# of admissions")

# load vector representation
dt.out <- read.table("old_output.txt", sep=" ", skip=2)
dt.out[,ncol(dt.out)] <- NULL
# matrix of icd codes with row names
row.names(dt.out) <- dt.out[,1]
dt.out <- dt.out[,-1]

# PCM: Primary code matching
# HDM: Hamming distance matching
# SM: Sequential matching
# WVM: Word2Vec matching
# matrix of sum of icd codes
sumCodes <- function(x)
{
  x <- x[!is.na(x)]
  colSums(dt.out[x,])
}
set.seed(999) # set seed 100
ptm <- proc.time()
iter <- 180 # iterations 150
n.iter <- 0 # real iterations
N <- 200 # sample episodes 200
valid <- vector(mode="numeric", length=iter)
acc.pc <- vector(mode="numeric", length=iter)
acc.hd <- vector(mode="numeric", length=iter)
acc.sm <- vector(mode="numeric", length=iter)
acc.wc <- vector(mode="numeric", length=iter)
acc.we <- vector(mode="numeric", length=iter)
acc.wm <- vector(mode="numeric", length=iter)
acc.ws <- vector(mode="numeric", length=iter)
ir.trial <- vector(mode="numeric", length=iter)
ir.pc.match <- vector(mode="numeric", length=iter)
ir.hd.match <- vector(mode="numeric", length=iter)
ir.sm.match <- vector(mode="numeric", length=iter)
ir.wc.match <- vector(mode="numeric", length=iter)
ir.we.match <- vector(mode="numeric", length=iter)
ir.wm.match <- vector(mode="numeric", length=iter)
ir.ws.match <- vector(mode="numeric", length=iter)
# sample hospital and admit year
hos_year <- df %>% distinct(hos, year) %>% sample_n(iter)
for(x in 1:iter) {
  # episodes with same hospital and admit year
  df.hy <- df %>% filter(hos==hos_year$hos[x],year==hos_year$year[x])
  # sample each patient 1 episode
  trial <- df.hy %>% group_by(pid) %>% sample_n(1) %>% as.data.frame()
  # sample N episodes from distinct patients
  if(N>nrow(trial)) {
    # cannot find enough N episodes in this iteration
    next
  }
  rnd <- sample(nrow(trial),N,replace=F)
  trial <- trial[rnd,]
  # generate control group
  control <- df.hy %>% filter(!(df.hy$pid %in% trial$pid))
  # matching
  e.readmission.pc <- vector(mode="numeric", length=N)
  e.readmission.hd <- vector(mode="numeric", length=N)
  e.readmission.sm <- vector(mode="numeric", length=N)
  e.readmission.wc <- vector(mode="numeric", length=N)
  e.readmission.we <- vector(mode="numeric", length=N)
  e.readmission.wm <- vector(mode="numeric", length=N)
  e.readmission.ws <- vector(mode="numeric", length=N)
  real.trial <- trial
  real.trial$real <- -1
  real.match.pc <- data.frame()
  real.match.hd <- data.frame()
  real.match.sm <- data.frame()
  real.match.wc <- data.frame()
  real.match.we <- data.frame()
  real.match.wm <- data.frame()
  real.match.ws <- data.frame()
  icds_trial <- trial[,icd_col:ncol(trial)]
  icds_trial$len <- rowSums(!is.na(icds_trial))
  for(i in 1:nrow(trial)) {
    # generate validation 
    validation <- control %>% 
      filter(control$age==as.character(trial[i,"age"]) & 
               control$gender==as.character(trial[i,"gender"]))
    if(nrow(validation)==0) {
      # cannot find validation group for this patient
      next
    }
    real.trial$real[i] <- i
    adm_tri <- trial$readm[i] %>% as.character()
    len_tri <- icds_trial$len[i]
    S_validation <- validation %>% select(sepdate,dth_date,readm,icd_col:ncol(validation))
    # 1. primary code matching
    # match the first ICD code
    S <- S_validation %>% filter(S_validation[,4]==as.character(icds_trial[i,1]))
    # can find matching cases
    if(nrow(S)>0) {
      mc <- S[sample(nrow(S),1),]
      readm <- mc[1,"readm"] %>% as.character() 
      e.readmission.pc[i] <- ifelse(adm_tri==readm,1,0)
      real.match.pc <- rbind(real.match.pc,mc)
    }
    # cannot find any matching case
    if(nrow(S)==0) {
      # select randomly from validation
      mc <- S_validation[sample(nrow(S_validation),1),]
      readm <- mc[1,"readm"] %>% as.character()
      e.readmission.pc[i] <- ifelse(adm_tri==readm,1,0)
      real.match.pc <- rbind(real.match.pc,mc)
    }
    # 2. hamming distance
    icds_validation <- validation[,icd_col:ncol(validation)]
    icds_validation$len <- rowSums(!is.na(icds_validation))
    hd <- sapply(1:nrow(S_validation), function(y) {
      len_y <- icds_validation$len[y]
      if(len_y<len_tri) {
        total <- 0
        for(j in 1:len_y) {
          if(as.character(icds_trial[i,j])==as.character(icds_validation[y,j])) {
            total <- total + 1
          }
        }
      }
      if(len_y>=len_tri) {
        total <- 0
        for(j in 1:len_tri) {
          if(as.character(icds_trial[i,j])==as.character(icds_validation[y,j])) {
            total <- total + 1
          }
        }
      }
      return(total)
    })
    mc <- S_validation[which(hd==max(hd)),]
    mc <- mc[sample(nrow(mc),1),]
    readm <- mc[1,"readm"] %>% as.character() 
    e.readmission.hd[i] <- ifelse(adm_tri==readm,1,0)
    real.match.hd <- rbind(real.match.hd,mc)
    # 3. sequentially matching & 4. word2vec matching
    S <- S_validation
    matched <- data.frame()
    for(j in 1:len_tri) {
      S <- S %>% filter(S[,j+3]==as.character(icds_trial[i,j]))
      if(nrow(S)>0) {
        matched <- S
      }
      else {
        break
      }
    }
    # exact matched
    if(j==len_tri & nrow(S)>0) {
      mc <- matched[sample(nrow(matched),1),]
      readm <- mc[1,"readm"] %>% as.character() 
      # SM and WVM are the same
      e.readmission.sm[i] <- ifelse(adm_tri==readm,1,0)
      e.readmission.wc[i] <- e.readmission.sm[i]
      e.readmission.we[i] <- e.readmission.sm[i]
      e.readmission.wm[i] <- e.readmission.sm[i]
      e.readmission.ws[i] <- e.readmission.sm[i]
      real.match.sm <- rbind(real.match.sm,mc)
      real.match.wc <- rbind(real.match.wc,mc)
      real.match.we <- rbind(real.match.we,mc)
      real.match.wm <- rbind(real.match.wm,mc)
      real.match.ws <- rbind(real.match.ws,mc)
    }
    # not exact matched
    if(j<=len_tri & nrow(S)==0) {
      # matched at least 1 icd code
      if(nrow(matched)>0) {
        # SM: random match
        mc <- matched[sample(nrow(matched),1),]
        readm <- mc[1,"readm"] %>% as.character()
        e.readmission.sm[i] <- ifelse(adm_tri==readm,1,0)
        real.match.sm <- rbind(real.match.sm,mc)
        # WVM: find similar icd in the matched icds
        v.icds_tri <- dt.out[as.character(icds_trial[i,j]),]
        tmp <- matched %>% select(j+3) %>% t(.) 
        if(is.na(tmp)) { # |icd codes| of trial ep. > |icd codes| of validation ep.
          e.readmission.wc[i] <- e.readmission.sm[i]
          e.readmission.we[i] <- e.readmission.sm[i]
          e.readmission.wm[i] <- e.readmission.sm[i]
          e.readmission.ws[i] <- e.readmission.sm[i]
          real.match.wc <- rbind(real.match.wc,mc)
          real.match.we <- rbind(real.match.we,mc)
          real.match.wm <- rbind(real.match.wm,mc)
          real.match.ws <- rbind(real.match.ws,mc)
        }
        else {
          tmp <- dt.out[tmp,]
          # convert to matrices
          v.icds_tri <- t(as.matrix(v.icds_tri))
          tmp <- as.matrix(tmp)
          # cosine similarity
          # compute unit vectors
          u.icds_tri <- sqrt(sum(v.icds_tri^2))
          u.tmp <- sqrt(rowSums(tmp^2))
          # compute cosine similarity matrix
          Dc <- tmp %*% v.icds_tri
          Dc <- Dc / u.tmp 
          Dc <- Dc / u.icds_tri
          # estimated readmission
          mc <- matched[which.max(Dc),]
          readm <- mc[1,"readm"] %>% as.character()
          e.readmission.wc[i] <- ifelse(adm_tri==readm,1,0)
          real.match.wc <- rbind(real.match.wc,mc)
          # euclidean distance
          v.icds_tri <- t(v.icds_tri)
          # compute square of unit vectors
          uu.icds_tri <- sum(v.icds_tri^2)
          uu.tmp <- rowSums(tmp^2)
          De <- matrix(rep(uu.icds_tri, nrow(tmp)), nrow=1) 
          De <- De + matrix(rep(uu.tmp, 1), nrow=1, byrow=TRUE)
          De <- De - 2 * tcrossprod(v.icds_tri,tmp)
          # estimated readmission
          mc <- matched[which.min(De),]
          readm <- mc[1,"readm"] %>% as.character()
          e.readmission.we[i] <- ifelse(adm_tri==readm,1,0)
          real.match.we <- rbind(real.match.we,mc)
          # manhattan distance
          v.icds_tri <- t(v.icds_tri)
          Dm <- matrix(rep(v.icds_tri, each=nrow(tmp)),ncol = length(v.icds_tri))
          Dm <- rowSums(abs(Dm - tmp)) %>% as.numeric()
          # estimated readmission
          mc <- matched[which.min(Dm),]
          readm <- mc[1,"readm"] %>% as.character()
          e.readmission.wm[i] <- ifelse(adm_tri==readm,1,0)
          real.match.wm <- rbind(real.match.wm,mc)
        }
      }
      # no matched for 1st icd code
      if(nrow(matched)==0) {
        # SM: select randomly from validation
        mc <- S_validation[sample(nrow(S_validation),1),]
        readm <- mc[1,"readm"] %>% as.character()
        e.readmission.sm[i] <- ifelse(adm_tri==readm,1,0)
        real.match.sm <- rbind(real.match.sm,mc)
        # WVM: find similar icd in the whole validation
        # get 1st icd
        v.icds_tri <- dt.out[as.character(icds_trial[i,1]),]
        tmp <- validation %>% select(icd_col) %>% t(.)
        tmp <- dt.out[tmp,]
        # convert to matrices
        v.icds_tri <- t(as.matrix(v.icds_tri))
        tmp <- as.matrix(tmp)
        # cosine similarity
        # compute unit vectors
        u.icds_tri <- sqrt(sum(v.icds_tri^2))
        u.tmp <- sqrt(rowSums(tmp^2))
        # compute cosine similarity matrix
        Dc <- tmp %*% v.icds_tri
        Dc <- Dc / u.tmp 
        Dc <- Dc / u.icds_tri
        # estimated readmission
        mc <- S_validation[which.max(Dc),]
        readm <- mc[1,"readm"] %>% as.character()
        e.readmission.wc[i] <- ifelse(adm_tri==readm,1,0)
        real.match.wc <- rbind(real.match.wc,mc)
        # euclidean distance
        v.icds_tri <- t(v.icds_tri)
        # compute square of unit vectors
        uu.icds_tri <- sum(v.icds_tri^2)
        uu.tmp <- rowSums(tmp^2)
        De <- matrix(rep(uu.icds_tri, nrow(tmp)), nrow=1) 
        De <- De + matrix(rep(uu.tmp, 1), nrow=1, byrow=TRUE)
        De <- De - 2 * tcrossprod(v.icds_tri,tmp)
        # estimated readmission
        mc <- S_validation[which.min(De),]
        readm <- mc[1,"readm"] %>% as.character()
        e.readmission.we[i] <- ifelse(adm_tri==readm,1,0)
        real.match.we <- rbind(real.match.we,mc)
        # mathattan distance
        v.icds_tri <- t(v.icds_tri)
        Dm <- matrix(rep(v.icds_tri, each=nrow(tmp)),ncol = length(v.icds_tri))
        Dm <- rowSums(abs(Dm - tmp)) %>% as.numeric()
        # estimated readmission
        mc <- S_validation[which.min(Dm),]
        readm <- mc[1,"readm"] %>% as.character()
        e.readmission.wm[i] <- ifelse(adm_tri==readm,1,0)
        real.match.wm <- rbind(real.match.wm,mc)
      }
    }
    # 5. sum matching
    v.icds_tri <- sumCodes(icds_trial[i,1:icds_trial$len[i]])
    tmp <- t(apply(validation[,c(icd_col:ncol(validation))], 1, FUN=sumCodes))
    # convert to matrices
    v.icds_tri <- as.matrix(v.icds_tri)
    tmp <- as.matrix(tmp)
    # cosine similarity
    # compute unit vectors
    u.icds_tri <- sqrt(sum(v.icds_tri^2))
    u.tmp <- sqrt(rowSums(tmp^2))
    # compute cosine similarity matrix
    Dc <- tmp %*% v.icds_tri
    Dc <- Dc / u.tmp
    Dc <- Dc / u.icds_tri
    # estimated readmission
    mc <- S_validation[which.max(Dc),]
    readm <- mc[1,"readm"] %>% as.character()
    e.readmission.ws[i] <- ifelse(adm_tri==readm,1,0)
    real.match.ws <- rbind(real.match.ws,mc)
  }
  n.iter <- n.iter + 1
  # compute accuracy
  real.trial <- real.trial %>% filter(real!=-1)
  valid[x] <- nrow(real.trial)
  acc.pc[x] <- sum(e.readmission.pc)/valid[x]
  acc.hd[x] <- sum(e.readmission.hd)/valid[x]
  acc.sm[x] <- sum(e.readmission.sm)/valid[x]
  acc.wc[x] <- sum(e.readmission.wc)/valid[x]
  acc.we[x] <- sum(e.readmission.we)/valid[x]
  acc.wm[x] <- sum(e.readmission.wm)/valid[x]
  acc.ws[x] <- sum(e.readmission.ws)/valid[x]
  # compute IR of trial
  n.event <- sum(!is.na(real.trial[,"dth_date"]))
  real.trial[,"dth_date"][is.na(real.trial[,"dth_date"])] <- cs_date
  ir.trial[x] <- n.event/sum(as.numeric(as.duration(real.trial[,"dth_date"]-
                                                      real.trial[,"sepdate"])/(3600*24*365)))
  # compute IR of PCM
  n.event <- sum(!is.na(real.match.pc[,"dth_date"]))
  real.match.pc[,"dth_date"][is.na(real.match.pc[,"dth_date"])] <- cs_date
  ir.pc.match[x] <- n.event/sum(as.numeric(as.duration(real.match.pc[,"dth_date"]-
                                                         real.match.pc[,"sepdate"])/(3600*24*365)))
  # compute IR of HDM
  n.event <- sum(!is.na(real.match.hd[,"dth_date"]))
  real.match.hd[,"dth_date"][is.na(real.match.hd[,"dth_date"])] <- cs_date
  ir.hd.match[x] <- n.event/sum(as.numeric(as.duration(real.match.hd[,"dth_date"]-
                                                         real.match.hd[,"sepdate"])/(3600*24*365)))
  # compute IR of SM
  n.event <- sum(!is.na(real.match.sm[,"dth_date"]))
  real.match.sm[,"dth_date"][is.na(real.match.sm[,"dth_date"])] <- cs_date
  ir.sm.match[x] <- n.event/sum(as.numeric(as.duration(real.match.sm[,"dth_date"]-
                                                         real.match.sm[,"sepdate"])/(3600*24*365)))
  # compute IR of WVM + cosine
  n.event <- sum(!is.na(real.match.wc[,"dth_date"]))
  real.match.wc[,"dth_date"][is.na(real.match.wc[,"dth_date"])] <- cs_date
  ir.wc.match[x] <- n.event/sum(as.numeric(as.duration(real.match.wc[,"dth_date"]-
                                                         real.match.wc[,"sepdate"])/(3600*24*365)))
  # compute IR of WVM + euclidean
  n.event <- sum(!is.na(real.match.we[,"dth_date"]))
  real.match.we[,"dth_date"][is.na(real.match.we[,"dth_date"])] <- cs_date
  ir.we.match[x] <- n.event/sum(as.numeric(as.duration(real.match.we[,"dth_date"]-
                                                         real.match.we[,"sepdate"])/(3600*24*365)))
  # compute IR of WVM + manhattan
  n.event <- sum(!is.na(real.match.wm[,"dth_date"]))
  real.match.wm[,"dth_date"][is.na(real.match.wm[,"dth_date"])] <- cs_date
  ir.wm.match[x] <- n.event/sum(as.numeric(as.duration(real.match.wm[,"dth_date"]-
                                                         real.match.wm[,"sepdate"])/(3600*24*365)))
  # compute IR of WVM + sum
  n.event <- sum(!is.na(real.match.ws[,"dth_date"]))
  real.match.ws[,"dth_date"][is.na(real.match.ws[,"dth_date"])] <- cs_date
  ir.ws.match[x] <- n.event/sum(as.numeric(as.duration(real.match.ws[,"dth_date"]-
                                                         real.match.ws[,"sepdate"])/(3600*24*365)))
}

# compute accuracy of PCM
sum(acc.pc)/n.iter
# compute standard error
sd(acc.pc[acc.pc>0])/sqrt(length(acc.pc[acc.pc>0]))
# compute IR error of PCM
err.pc <- abs(ir.trial-ir.pc.match)
sum(err.pc)/n.iter
# compute standard error
sd(err.pc[err.pc>0])/sqrt(length(err.pc[err.pc>0]))

# compute accuracy of HDM
sum(acc.hd)/n.iter
# compute standard error
sd(acc.hd[acc.hd>0])/sqrt(length(acc.hd[acc.hd>0]))
# compute IR error of HDM
err.hd <- abs(ir.trial-ir.hd.match)
sum(err.hd)/n.iter
# compute standard error
sd(err.hd[err.hd>0])/sqrt(length(err.hd[err.hd>0]))

# compute accuracy of SM
sum(acc.sm)/n.iter
# compute standard error
sd(acc.sm[acc.sm>0])/sqrt(length(acc.sm[acc.sm>0]))
# compute IR error of SM
err.sm <- abs(ir.trial-ir.sm.match)
sum(err.sm)/n.iter
# compute standard error
sd(err.sm[err.sm>0])/sqrt(length(err.sm[err.sm>0]))

# compute accuracy of WVM + cosine
sum(acc.wc)/n.iter
# compute standard error
sd(acc.wc[acc.wc>0])/sqrt(length(acc.wc[acc.wc>0]))
# compute IR error of WVM + cosine
err.wc <- abs(ir.trial-ir.wc.match)
sum(err.wc)/n.iter
# compute standard error
sd(err.wc[err.wc>0])/sqrt(length(err.wc[err.wc>0]))

# compute accuracy of WVM + euclidean
sum(acc.we)/n.iter
# compute standard error
sd(acc.we[acc.we>0])/sqrt(length(acc.we[acc.we>0]))
# compute IR error of WVM + euclidean
err.we <- abs(ir.trial-ir.we.match)
sum(err.we)/n.iter
# compute standard error
sd(err.we[err.we>0])/sqrt(length(err.we[err.we>0]))

# compute accuracy of WVM + manhattan
sum(acc.wm)/n.iter
# compute standard error
sd(acc.wm[acc.wm>0])/sqrt(length(acc.wm[acc.wm>0]))
# compute IR error of WVM + manhattan
err.wm <- abs(ir.trial-ir.wm.match)
sum(err.wm)/n.iter
# compute standard error
sd(err.wm[err.wm>0])/sqrt(length(err.wm[err.wm>0]))

# compute accuracy of WVM + sum
sum(acc.ws)/n.iter
# compute standard error
sd(acc.ws[acc.ws>0])/sqrt(length(acc.ws[acc.ws>0]))
# compute IR error of WVM + sum
err.ws <- abs(ir.trial-ir.ws.match)
sum(err.ws)/n.iter
# compute standard error
sd(err.ws[err.ws>0])/sqrt(length(err.ws[err.ws>0]))
proc.time() - ptm


library(ggplot2)
n.run <- 150
ir_stats <- cbind(c(1:n.run),rep("0. Case (ground truth)",n.run),sort(ir.trial[ir.trial>0][1:n.run],
                                                             decreasing=T)) %>%
  rbind(.,cbind(c(1:n.run),rep("1. PCM",n.run),sort(ir.pc.match[ir.pc.match>0][1:n.run],
                                                         decreasing=T))) %>%
  rbind(.,cbind(c(1:n.run),rep("2. HDM",n.run),sort(ir.hd.match[ir.hd.match>0][1:n.run],
                                                          decreasing=T))) %>%
  #rbind(.,cbind(c(1:n.run),rep("3. SM",n.run),sort(ir.sm.match[ir.sm.match>0][1:n.run],
  #                                                      decreasing=T))) %>%
  rbind(.,cbind(c(1:n.run),rep("3. CSM",n.run),sort(ir.ws.match[ir.ws.match>0][1:n.run],
                                                                decreasing=T))) %>%
  rbind(.,cbind(c(1:n.run),rep("4. WVM + cosine",n.run),sort(ir.wc.match[ir.wc.match>0][1:n.run],
                                                                  decreasing=T))) %>%
  rbind(.,cbind(c(1:n.run),rep("5. WVM + euclidean",n.run),sort(ir.we.match[ir.we.match>0][1:n.run],
                                                                     decreasing=T))) %>%
  rbind(.,cbind(c(1:n.run),rep("6. WVM + manhattan",n.run),sort(ir.wm.match[ir.wm.match>0][1:n.run],
                                                                     decreasing=T))) %>%
  as.data.frame()
colnames(ir_stats) <- c("iter","Incidence_rate_of","IR")
ir_stats$iter <- as.numeric(as.character(ir_stats$iter))
ir_stats$Incidence_rate_of <- as.character(ir_stats$Incidence_rate_of)
ir_stats$IR <- as.numeric(as.character(ir_stats$IR))

blank_theme <- theme_bw() +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.grid=element_blank(),
    axis.ticks.length = unit(.5, "cm"),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.key.size = unit(0.8, "cm"),
    legend.text = element_text(size = 20),
    legend.justification = c(1, 1), 
    legend.position = c(1, 1)
  )

ir_stats %>% ggplot(aes(x = iter, y = IR, color = Incidence_rate_of)) + 
  geom_line(aes(group=Incidence_rate_of), size = 1) +
  xlab("Iteration") + ylab("Incidence rate") + blank_theme
  #ggtitle("Incidence Rate Comparison") +
  #theme(legend.justification = c(1, 1), legend.position = c(1, 1))



