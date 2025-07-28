# demo.R - Demo script for MP running
# jmMSE/demo.R

# Copyright (c) WMR, 2025.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


# INSTALL only if needed
# install.packages(c("TAF", "qs", "FLjjm", "FLCore", "FLFishery",
#   "FLasher", "mse", "ggplotFL", "mseviz"),
#  repos = c('https://sprfmo.r-universe.dev', 'https://cloud.r-project.org'))

# SETUP

source(here::here("R","config.R") )
mkdir("demo")

# LOAD H1 inputs: om, oem, iem
qload('data/h1_1.07.qs')

# SET OEM
method(oem) <- shortcut.oem
args(oem)$jjms <- FALSE

# SELECT stock for MP
mseargs$stock <- 1
# plan(multisession, workers=2)
plan(sequential)
#plan(multicore, workers=5)

# --- good SHORTCUT (model_shortcut_buffer.R) {{{ ---

# metric devs, CV=15%, rho=0
metdevs <- rlnormar1(dims(om)$iter, sdlog=0.15, rho=0, years=seq(1970, 2045))

# MP control
ctrl <- mpCtrl(list(
  # EST
  est = mseCtrl(method=shortcut.sa, args=list(metric="depletion",
                                              devs=metdevs, B0=refpts(om)$SB0)),
  # HCR
  hcr = mseCtrl(method=buffer.hcr,
                args=list(target=1000, bufflow=0.30, buffupp=0.50, lim=0.10, min=0,
                          metric="depletion")),
  # ISYS
  isys=mseCtrl(method=split.is, args=list(split=catch_props(om)$last5))
))

run <- mp(om, oem, iem=iem, ctrl=ctrl, args=mseargs)

# RUN2
# MP control
ctrl <- mpCtrl(list(
  # EST
  est = mseCtrl(method=shortcut.sa, args=list(metric="depletion",
                                              devs=metdevs, B0=refpts(om)$SB0)),
  # HCR
  hcr = mseCtrl(method=buffer.hcr,
                args=list(target=1000, bufflow=0.20, buffupp=0.20, lim=0.10, min=0,
                          metric="depletion")),
  # ISYS
  isys=mseCtrl(method=split.is, args=list(split=catch_props(om)$last5))
))
# RUN2
run2 <- mp(om, oem, iem=iem, ctrl=ctrl, args=mseargs)

#---RUN3 ----
# MP control
ctrl <- mpCtrl(list(
  # EST
  est = mseCtrl(method=shortcut.sa, args=list(metric="depletion",
                                              devs=metdevs, B0=refpts(om)$SB0)),
  # HCR
  hcr = mseCtrl(method=buffer.hcr,
                args=list(target=1400, bufflow=0.20, buffupp=0.20, lim=0.10, min=0,
                          metric="depletion")),
  # ISYS
  isys=mseCtrl(method=split.is, args=list(split=catch_props(om)$last5))
))
# RUN3 
run3 <- mp(om, oem, iem=iem, ctrl=ctrl, args=mseargs)
# PLOT
plot(om, run)

source(here::here("R","extract_tracking.R"))
# Example usage for all metrics:
df_all <- extract_tracking_array(run)
df_all$run<-"Run1"
df_tmp <- extract_tracking_array(run2)
df_tmp$run<-"Run2"
df_tm3 <- extract_tracking_array(run3)
df_tm3$run<-"Run3"
 summary(df_all)
 names(df_all)

df_all <- rbind(df_all,df_tmp, df_tm3)
# # Example usage:
df_all |> dplyr::filter(SB.om<.50e5) |> ggplot(aes(x=(SB.om),y=(tac.hcr), color=run)) + 
geom_point(size=1.0,alpha=.5) + ggthemes::theme_few() + geom_smooth()

runs <- FLmses(list(run1=run, run2=run2, run3=run3 ), statistics=statistics,
               years=2024:2045, metrics=mets, type="scgood_buffer")

writePerformance(performance(runs), file="demo/runs.dat.gz")
performance(runs)[statistic == 'green' & year %in% 2034:2042, .(green=mean(data))]

perf  <- readPerformance("demo/runs.dat.gz")
omperf <- performance(om, statistics=statistics[c("C","F","SB")])
names(perf)
names(omperf)
plot_omperf(omperf)

plot_mp_panels(perf)
plot_C_SSB_F_stack(perf)
plotSummary(perf)
plot_profile_targets_stack(perf)
plot_profile_metric(perf)


#--Now try a cpue again
method(oem) <- sampling.oem

imu <- yearMeans(observations(oem)$CJM$idx[[3]]@index)
isd <- sqrt(yearVars(observations(oem)$CJM$idx[[3]]@index))
isd/imu

ctrl <- mpCtrl(list(
  # EST
  est = mseCtrl(method=cpuescore2.ind, args=list(index=3, mean=imu, sd=isd)),
  # HCR
  hcr = mseCtrl(method=buffer.hcr,
    args=list(target=1000, bufflow=-1, buffupp=1, lim=-2, min=0,
      sloperatio=0.15, metric="zscore")),
  # ISYS
  isys=mseCtrl(method=split.is, args=list(split=catch_props(om)$last5))
))

# TEST:
tes <- mp(om, oem, iem=iem, ctrl=ctrl, args=mseargs)
plot(om,tes)
df_tes <- extract_tracking_array(tes)  |> dplyr::mutate(run="CPUE3")
names(df_tes)
df_tes |> dplyr::filter(SB.om<.50e5) |> ggplot(aes(x=(SB.om),y=(tac.hcr), color=run)) + 
geom_point(size=1.0,alpha=.5) + ggthemes::theme_few() + geom_smooth()

df_tes |> dplyr::filter(SB.om<.50e5) |> ggplot(aes(x=(cpue.ind),y=(tac.hcr), color=run)) + 
geom_point(size=1.0,alpha=.5) + ggthemes::theme_few() + geom_smooth()

#---tes2----------
ctrl <- mpCtrl(list(
  # EST
  est = mseCtrl(method=cpuescore2.ind, args=list(index=3, mean=imu, sd=isd)),
  # HCR
  hcr = mseCtrl(method=buffer.hcr,
    args=list(target=1000, bufflow=-.5, buffupp=.5, lim=-2, min=0,
      sloperatio=0.15, metric="zscore")),
  # ISYS
  isys=mseCtrl(method=split.is, args=list(split=catch_props(om)$last5))
))

# TEST:
tes2 <- mp(om, oem, iem=iem, ctrl=ctrl, args=mseargs)

#---tes3----------
ctrl <- mpCtrl(list(
  # EST
  est = mseCtrl(method=cpuescore2.ind, args=list(index=3, mean=imu, sd=isd)),
  # HCR
  hcr = mseCtrl(method=buffer.hcr,
    args=list(target=1400, bufflow=-.5, buffupp=.5, lim=-2, min=0,
      sloperatio=0.15, metric="zscore")),
  # ISYS
  isys=mseCtrl(method=split.is, args=list(split=catch_props(om)$last5))
))

# TEST:
tes3 <- mp(om, oem, iem=iem, ctrl=ctrl, args=mseargs)

#---tes4----------
ctrl <- mpCtrl(list(
  # EST
  est = mseCtrl(method=cpuescorejim.ind, args=list(index=3, mean=imu, sd=isd)),
  # HCR
  hcr = mseCtrl(method=buffer.hcr,
    args=list(target=1400, bufflow=-.5, buffupp=.5, lim=-2, min=0,
      sloperatio=0.15, metric="zscore")),
  # ISYS
  isys=mseCtrl(method=split.is, args=list(split=catch_props(om)$last5))
))

# TEST:
tes4 <- mp(om, oem, iem=iem, ctrl=ctrl, args=mseargs)

df_tes <- rbind(
extract_tracking_array(tes)  |> dplyr::mutate(run="i3_1.0 1.0"),
extract_tracking_array(tes2)  |> dplyr::mutate(run="i3_0.5 1.0"),
extract_tracking_array(tes3)  |> dplyr::mutate(run="i3_0.5 1.4"),
extract_tracking_array(tes4)  |> dplyr::mutate(run="i3_0.5 Jim")
)
names(df_tes)
df_tes |> dplyr::filter(SB.om<.50e5) |> ggplot(aes(x=(SB.om),y=(tac.hcr), color=run)) + 
geom_point(size=1.0,alpha=.5) + ggthemes::theme_few() + geom_smooth()

df_tes |> dplyr::filter(SB.om<.30e5) |> ggplot(aes(x=(cpue.ind),y=(tac.hcr), color=run)) + 
geom_point(size=1.0,alpha=.5) + ggthemes::theme_few() + geom_smooth()

df_tes |> dplyr::filter(SB.om<.30e5) |> ggplot(aes(x=(cpue.ind),y=(tac.hcr), color=run)) + 
geom_point(size=1.0,alpha=.5) + ggthemes::theme_few() + geom_smooth()

names(df_tes)



# TODO: labels
labels <- list(
  # OMs
  h1_1.07_='H1', h1_1.01_enso_='H1EN',
  # est
  scgood='SG', scmedium='SM', scbad='SB',
  cpue2='C2', cpue3='C3', cpue36='C36',
  # hcr
  buffer='',
  # run
  hcr_target='target', hcr_bufflow='bufflow', tune='K',
  # clean up
  `__`='_')

for(i in seq(labels)) {
  perf[, mp:=gsub(names(labels)[i], labels[i], mp)]
}

# - Slick
sli <- getSlick(perf, omperf, kobeyrs=2034:2042)

# RUN app



# - EXPLORE HCR targets

targets <- mps(om, oem, ctrl=ctrl, args=mseargs,
               hcr=list(target=seq(500, 1500, length=3)))

# PLOT
plot(om, targets)

# COMPUTE & ASSIGN performance
performance(targets) <- performance(targets, statistics=statistics, 
                                    type="scgood_buffer")

# STORE performance table
writePerformance(performance(targets), file="demo/performance.dat.gz", overwrite=TRUE)
readPerformance(performance(targets), file="demo/performance.dat.gz", overwrite=TRUE)

# SAVE runs
qsave(targets, file="demo/h1_scgood_buffer_targets.qs")

# - EXPLORE HCR low buffer
bufflows <- mps(om, oem, ctrl=ctrl, args=mseargs,
                hcr=list(bufflow=seq(0.2, 0.40, length=3)))

# PLOT
plot(om, bufflows)

# COMPUTE & ASSIGN performance
performance(bufflows) <- performance(bufflows, statistics=statistics,
                                     type="scgood_buffer")

# STORE performance table
writePerformance(performance(bufflows), file="demo/performance.dat.gz")

# SAVE runs
qsave(bufflows, file="demo/h1_scgood_buffer_bufflows.qs")
bufflows <- qread(file="demo/h1_scgood_buffer_bufflows.qs")

# - tune for P(Kobe=green) = 60%

system.time(
  tune06 <- tunebisect(om, oem=oem, iem=iem, control=ctrl, args=mseargs, 
                       statistic=statistics["green"], years=2034:2042,
                       tune=list(target=c(500, 2000)), prob=0.6, tol=0.01, maxit=12)
)

tune <- FLmses(list(tune06=tune06), statistics=statistics,
               years=2024:2045, metrics=mets, type="scgood_buffer")

# STORE performance tabl
writePerformance(performance(tune), file="demo/performancetune.dat.gz")

# SAVE run
qsave(tune, file="demo/h1_scgood_buffer_tune.qs")

# INSPECT results

plot(om, tune)
# Has tuning probability been achieved?
performance(tune)[statistic == 'green' & year %in% 2034:2042, .(green=mean(data))]

# How is P(Kobe=green) changing in time?
performance(tune)[statistic == 'green', .(green=mean(data)), by=year]

# }}}

# --- SLICK (output.R) {{{ ----

library(Slick)

perf <- readPerformance("demo/performance.dat.gz")
omperf <- performance(om, statistics=statistics[c("C","F","SB")])

# TODO: labels
labels <- list(
  # OMs
  h1_1.07_='H1', h1_1.01_enso_='H1EN',
  # est
  scgood='SG', scmedium='SM', scbad='SB',
  cpue2='C2', cpue3='C3', cpue36='C36',
  # hcr
  buffer='',
  # run
  hcr_target='target', hcr_bufflow='bufflow', tune='K',
  # clean up
  `__`='_')

for(i in seq(labels)) {
  perf[, mp:=gsub(names(labels)[i], labels[i], mp)]
}

# - Slick
sli <- getSlick(perf, omperf, kobeyrs=2034:2042)

# RUN app
App(slick=sli)

# }}}

# --- PLOT results (report.R) {{{ ----

library(mseviz)

# COMPUTE averages across periods

periods <- list(tuning=2034:2045, short=2025:2027, medium=2027:2032, long=2033:2045)

perf <- periodsPerformance(perf, periods)

# PLOT TS

# PLOT BPs  TODO: ylim
plotBPs(perf[period == 'tuning'])

# PLOT TOs
plotTOs(perf[period == 'tuning'])

# - TABLES (short, medium, long)

library(flextable)

as_flextable(summTable(perf[period == "tuning"]))

# }}}

