# utilities.R - Extra functions
# jmMSE/utilities.R

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


# controlAssign {{{

controlAssign <- function(target, years, quant = "fbar", ratio) {
  if (quant %in% c("f", "fbar")) {
    trg <- lapply(years, function(y) {
      list(
        year = y, biol = 1, minAge = 1, maxAge = 12,
        quant = quant, value = c(target)
      )
    })
  } else {
    trg <- lapply(years, function(y) {
      list(
        year = y, biol = 1,
        quant = quant, value = unlist(target)
      )
    })
  }

  fwdControl(c(
    trg,
    lapply(years, function(y) {
      list(
        year = y, relYear = y, fishery = 1, catch = 1,
        relFishery = 2, relCatch = 1, quant = "catch", value = ratio[1]
      )
    }),
    lapply(years, function(y) {
      list(
        year = y, relYear = y, fishery = 3, catch = 1,
        relFishery = 2, relCatch = 1, quant = "catch", value = ratio[3]
      )
    }),
    lapply(years, function(y) {
      list(
        year = y, relYear = y, fishery = 4, catch = 1,
        relFishery = 2, relCatch = 1, quant = "catch", value = ratio[4]
      )
    })
  ))
}
# }}}

# fixidxcjm {{{

fixidxcjm <- function(oem) {
  #  SHORT old indices TODO: MOVE to fwdWindow(FLoem)

  if ("CJM" %in% names(observations(oem))) {
    observations(oem)$CJM$idx$Chile_AcousCS <-
      window(observations(oem)$CJM$idx$Chile_AcousCS, end = 2009)

    observations(oem)$CJM$idx$DEPM <-
      window(observations(oem)$CJM$idx$DEPM, end = 2008)

    observations(oem)$CJM$idx$Peru_Acoustic <-
      window(observations(oem)$CJM$idx$Peru_Acoustic, end = 2013)

    #  SET sel.pattern to last 20 years
    for (i in c("Chile_AcousN", "Chile_CPUE", "Peru_CPUE", "Offshore_CPUE")) {
      index.q(observations(oem)$CJM$idx[[i]])[, ac(2025:2045)] <-
        yearMeans(window(index.q(observations(oem)$CJM$idx[[i]]),
          start = 2001, end = 2020
        ))
      sel.pattern(observations(oem)$CJM$idx[[i]])[, ac(2025:2045)] <-
        yearMeans(window(sel.pattern(observations(oem)$CJM$idx[[i]]),
          start = 2001, end = 2020
        ))
    }
  } else {
    observations(oem)$Southern$idx$Chile_AcousCS <-
      window(observations(oem)$Southern$idx$Chile_AcousCS, end = 2009)

    observations(oem)$Southern$idx$DEPM <-
      window(observations(oem)$Southern$idx$DEPM, end = 2008)

    observations(oem)$North$idx$Peru_Acoustic <-
      window(observations(oem)$North$idx$Peru_Acoustic, end = 2013)
  }
  return(oem)
}
# }}}

# catch_props {{{

catch_props <- function(x) {
  catch_byfishery <- iter(window(catch(x, "fishery") / catch(x), end = 2024), 1)

  return(lapply(setNames(c(-5, -10, -15, -20),
    nm = paste0("last", c(5, 10, 15, 20))
  ), function(y) {
    unlist(lapply(window(catch_byfishery, start = y), yearMeans))
  }))
}
# }}}


# cjm.iem {{{
cjm.iem <- function(ctrl, args, tracking, F3prop, correction = NULL) {
  # VARIABLES
  my <- ac(args$ay + args$management_lag)

  # CORRECT total catch
  if (!is.null(correction)) {
    iters(ctrl)[, "value", ] <- iters(ctrl)[, "value", ] * (1 + correction)
  }

  # COMPUTE current F3 ratio
  ratio3 <- c(iters(ctrl)[3, 2, , drop = FALSE]) /
    c(apply(iters(ctrl)[, 2, , drop = FALSE], 2:3, sum))

  # MIMIC separate decision for fishery 3
  iters(ctrl)[3, 2, ] <- (c(F3prop[, my]) * iters(ctrl)[3, 2, ]) / ratio3

  # TRACK new total catch
  track(tracking, "tac.iem", my) <- c(apply(iters(ctrl)[, 2, , drop = FALSE], 2, sum))

  # TODO: TRACK catch by fishery

  return(list(ctrl = ctrl, tracking = tracking))
}

# }}}

# getSlick {{{

#' Convert mse output to a slick object
#'
#' @param df A data.table containing performance indicator values for each operating
#' model (OM), management procedure (MP) and simulation (sim)
#' @param kobey Range of years to use in Kobe plot, e.g. tuning period
#'
#' @return
#' @export
#'
#' @examples
getSlick <- function(x, y, kobeyrs = unique(x$year)) {
  # GET dims
  mps <- unique(x$mp)
  oms <- unique(x$om)

  minyr <- min(x$year)

  # CHECK year is numeric in both x & y
  x[, year := as.numeric(year)]
  y[, year := as.numeric(year)]

  # MULTIPLY y by MP
  ys <- rbindlist(lapply(setNames(nm = mps), function(i) {
    #  EXTRACT type & mp
    dat <- x[mp == i, .(type, run, mp)][1, ]
    # BIND to OM perf
    return(cbind(y[year < minyr], dat))
  }))

  # COMBINE MP and OM perfs
  tab <- rbindlist(list(ys, x), use.names = TRUE)

  # CHANGE columns class
  tab[, iter := as.numeric(iter)]
  tab[, year := as.numeric(year)]

  # SET order
  setorder(tab, statistic, mp, om, iter)

  # BUG: DROP large F (>4), FMSY (>8), IACC (>200)
  tab[statistic == "IACC" & data > 400, data := 400]
  tab[statistic == "F" & data > 4, data := 4]
  tab[statistic == "FMSY" & data > 8, data := 8]

  # CREATE Slick
  sli <- FLslick(tab, kobey = as.numeric(kobeyrs))

  # SET current year, end of OM
  sli@Timeseries@TimeNow <- as.numeric(minyr) # as.numeric(y[, max(year)])

  # CHANGES to Slick

  # - TimeSeries
  #   - [X] SOLVE year gap OM - MPs
  #   - [ ] ADD target and limit (F, C, SB)
  #   - FILTER statistics from OM period

  # - Trade Offs
  #   - or DROP outliers

  # - Kobe
  #   - Limits to show green
  #   - No green in time
  #   - Y axis to 200%
  sli@Kobe@Limit <- NULL
  sli@Kobe@Time <- as.numeric(kobeyrs)
  tmp <- sli@Kobe@Value
  sli@Kobe@Value[, , , 1, ] <- tmp[, , , 2, ]
  sli@Kobe@Value[, , , 2, ] <- tmp[, , , 1, ]
  sli@Kobe@Code <- c("SB/SBMSY", "FB/FBMSY")
  #   - Year range, not just start

  return(sli)
}


FLslick <- function(df, kobey = NULL) {
  df <- as.data.frame(df)

  # Basic checks ------------------------------------------------------------

  # Check if exists/empty/expected class
  if (is.null(df) || !is.data.frame(df) || nrow(df) <= 0 || !exists("df")) {
    stop("The input object is not valid")
  }

  # TODO: add check with names mapping?
  if (!class(df$iter) %in% "numeric") {
    df$iter <- as.numeric(df$iter)
  }


  # MPs ---------------------------------------------------------------------

  # tuneMP <-   unique(df$mp)[grep("tune", unique(df$mp), ignore.case = TRUE)]

  # df <- subset(df, mp %in% tuneMP)

  mps <- MPs(
    Code = unique(df$mp),
    # MATCH
    Label = unique(df$mp),
    Description = unique(df$mp)
  )

  # OMs ---------------------------------------------------------------------

  if (sum(names(df) == "om") > 0) {
    df$omDesc <- ifelse(grepl("h1", df$om), "Single stock",
      ifelse(grepl("h2", df$om), "Two stocks",
        ifelse(grepl("(?=.*h2)(?=.*m)", df$om, perl = TRUE), "Two stocks with movement rates", NA)
      )
    )

    # S stands for stock
    fac <- data.frame(
      Factor = rep("S", length(unique(df$om))),
      Level = unique(df$om),
      Description = unique(df$omDesc)
    )

    des <- data.frame(S = unique(df$om))
    oms <- OMs(
      Factors = fac,
      Design = des
    )
  } else {
    stop("Object does not have an om specified") # too strict? Warning and dummy col maybe?
  }

  # Determine cutoff year for filtering and timeseries ----------------------
  # FIND first year of MPs
  cutoff <- as.numeric(format(Sys.Date(), "%Y")) - 1


  # Keep different objects for respective plots -----------------------------

  # For all plots except the timeseries
  perf <- subset(df, df$year > cutoff)

  ## For kobe plot
  if (missing(kobey)) {
    warning("You have not specified the tuning years") # TODO better warning
    sbf <- subset(perf, statistic %in% c("SBMSY", "FMSY"))
  } else {
    sbf <- subset(perf, statistic %in% c("SBMSY", "FMSY") & year %in% kobey)
  }

  # For timeseries

  tim <- subset(df, statistic %in% c("C", "F", "SB"))

  # Boxplot -----------------------------------------------------------------
  # Dimensions: c(nsim, nOM, nMP, and nPI)
  perf <- within(perf, {
    iter <- ave(mp, list(mp, statistic, om), FUN = seq_along)
  })
  perf <- perf[order(perf$statistic, perf$mp, perf$om, perf$iter), ]


  dff <- array(
    data = perf$data,
    dim = c(
      length(unique(perf$iter)),
      length(unique(perf$om)),
      length(unique(perf$mp)),
      length(unique(perf$statistic))
    )
  )


  boxpl <- Boxplot(
    Code = unique(perf$statistic),
    Label = unique(perf$name),
    Description = unique(perf$desc),
    Value = dff,
    Preset = list(),
    Defaults = list("overall", "boxplot")
  )

  # slick <- Slick(MPs = mps)
  # Boxplot(slick) <- boxpl

  # plotBoxplot(slick)
  #
  #
  # # Kobe --------------------------------------------------------------------
  # # Dimensions: c(nsim, nOM, nMP, nPI, nTS)


  sbf <- sbf[order(sbf$year, sbf$name, sbf$mp, sbf$om, sbf$iter), ]

  dff <- array(
    data = sbf$data,
    dim = c(
      length(unique(sbf$iter)),
      length(unique(sbf$om)),
      length(unique(sbf$mp)),
      length(unique(sbf$name)),
      length(unique(sbf$year))
    )
  )
  # TODO change Time, needs to be 2034 to 2042
  kobe <- Kobe(
    Code = unique(sbf$statistic),
    Label = unique(sbf$name),
    Description = unique(sbf$desc),
    Time = numeric(length(unique(sbf$year))),
    TimeLab = "Year",
    Value = dff,
    Preset = list(),
    Target = 1,
    Limit = 1
  )

  #
  # # Kobe(slick) <- kobe
  #
  # # plotKobe(slick)
  #
  # # Quilt -------------------------------------------------------------------
  # # Dimensions: c(nOM, nMP, and nPI)
  #
  #
  # # Check if statistic is empty. If it is remove
  #
  perf <- perf[!apply(perf["statistic"], 1, function(x) all(is.na(x) | x == "")), ]

  # Mean or median?

  perfav <- setNames(aggregate(perf$data, by = list(perf$mp, perf$biol, perf$statistic, perf$name, perf$desc, perf$om), FUN = mean, na.rm = TRUE), c("mp", "biol", "statistic", "name", "desc", "om", "data"))


  perfav <- perfav[order(perfav$name, perfav$mp, perfav$om), ]


  dff <- array(
    data = perfav$data,
    dim = c(
      length(unique(perfav$om)),
      length(unique(perfav$mp)),
      length(unique(perfav$name))
    )
  )
  quilt <- Quilt(
    Code = unique(perfav$statistic),
    Label = unique(perfav$name),
    Description = unique(perfav$desc),
    Value = dff,
    Preset = list(),
    Color = c("white", "darkblue"),
    MinValue = 0,
    MaxValue = as.numeric(NA)
  )

  # # Set min and max values for the color coding
  # # MinValue(quilt) <- c(0,0,0,NA)
  # # MaxValue(quilt) <- c(1,1,1,NA)
  #
  # # slick <- Slick(MPs = mps)
  # # Quilt(slick) <- quilt
  # #
  # # plotQuilt(slick)
  #
  #
  #
  # # Spider ------------------------------------------------------------------
  # # Dimensions: c(nOM, nMP, and nPI)
  #
  # #scale values for non prob pi between 0 and 1
  #
  # # scale01 <- function(x) {
  # #   (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  # #

  # perfsp <- subset(perf, statistic %in% c( "green", "PSBMSY", "PSBlim", "PC0" ))
  perfsp <- perf[order(perf$name, perf$mp, perf$om), ]

  perfsp <- as.data.frame(data.table(perfsp)[, data := data / max(data), by = .(om, mp, statistic, name)])

  dffs <- array(
    data = perfsp$data,
    dim = c(
      length(unique(perfsp$om)),
      length(unique(perfsp$mp)),
      length(unique(perfsp$name))
    )
  )

  spider <- Spider(
    Code = unique(perfsp$statistic),
    Label = unique(perfsp$name),
    Description = unique(perfsp$desc),
    Value = dffs,
    Preset = list()
  )

  # # Check(slick)
  # #
  # # Spider(slick) <- spider
  #
  # # plotSpider(slick)
  #
  # # # Timeseries --------------------------------------------------------------
  # # # Dimensions: c(nsim, nOM, nMP, nPI, nTS)



  tim <- tim[order(tim$year, tim$name, tim$mp, tim$om, tim$iter), ]

  dff <- array(
    data = tim$data,
    dim = c(
      length(unique(tim$iter)),
      length(unique(tim$om)),
      length(unique(tim$mp)),
      length(unique(tim$name)),
      length(unique(tim$year))
    )
  )

  timeseries <- Timeseries(
    Code = unique(tim$statistic),
    Label = unique(tim$name),
    Description = unique(tim$desc),
    Time = c(
      rev(seq(cutoff, by = -1, length.out = cutoff - min(as.numeric(df$year)) + 1)),
      seq(cutoff + 1, by = 1, length.out = max(as.numeric(df$year)) - (cutoff))
    ),
    TimeNow = cutoff,
    TimeLab = "Year",
    Value = dff,
    Preset = list(),
    Target = NULL,
    Limit = NULL
  )
  #
  #
  # Check(timeseries)
  # slick <- Slick()
  # MPs(slick) <- mps
  # Timeseries(slick) <- timeseries
  #
  # plotTimeseries(slick)
  # #
  #
  #
  # # # Tradeoff ----------------------------------------------------------------
  # # # Dimensions: c(nOM, nMP, nPI)

  tr <- setNames(aggregate(perf$data, by = list(perf$mp, perf$statistic, perf$om), FUN = mean, na.rm = TRUE), c("mp", "statistic", "om", "data"))
  tr <- tr[order(tr$statistic, tr$mp, tr$om), ]


  dffs <- array(
    data = tr$data,
    dim = c(
      length(unique(tr$om)),
      length(unique(tr$mp)),
      length(unique(tr$statistic))
    )
  )

  tradeoff <- Tradeoff(
    Code = unique(tr$statistic),
    Label = unique(tr$statistic),
    Description = unique(tr$statistic),
    Value = dffs,
    Preset = list()
  )


  # Slick object ------------------------------------------------------------



  sl <- Slick(
    OMs = oms,
    MPs = mps,
    Boxplot = boxpl,
    Kobe = kobe,
    Quilt = quilt,
    Spider = spider,
    Timeseries = timeseries,
    Tradeoff = tradeoff
  )

  # Metadata intro ----------------------------------------------------------

  sl@Title <- "SPRFMO CJM MSE"
  sl@Date <- "February 2025"
  sl@Email <- "iagomosqueira@wur.nl"
  sl@Introduction <- "MSE analysis of candidate MPs for SPRFMO CJM"
  sl@Subtitle <- "COMM13 MSE workshop"
  sl@Author <- "Iago Mosqueira, Karolina Molla-Gazi, Justin Tiano"
  sl@Institution <- "WMR, The Netherlands"

  return(sl)
} # }}}

# OLD bank_borrow.is {{{

OLDbank_borrow.is <- function(
    stk, ctrl, args, split, bank = NULL, borrow = NULL,
    tracking) {
  # DIMS
  ay <- args$ay
  iy <- args$iy
  pys <- ac(seq(
    args$ay + args$management_lag,
    args$ay + args$management_lag + args$frq - 1
  ))

  # CHECK either bank or borrow, not both

  # APPLY rules for banking and/or borrowing
  # - stock status (or TAC changes)

  # DO bank
  if (!is.null(bank)) {
    if (is(bank, "FLQuant")) {
      bank <- c(propagate(bank[, pys], args$it))
    }

    # COMPUTE amount
    banking <- ctrl$value * bank

    # CORRECT current advice
    ctrl$value <- ctrl$value - banking

    # STORE banking value in decision year
    mse::track(tracking, "banking.isys", ay) <- banking

    # APPLY banking
    banked <- tracking[[1]]["banking.isys", ac(ay - 1)]

    if (!all(is.na(banked))) {
      ctrl$value <- ctrl$value + banked
    }
  }

  # DO borrow
  if (!is.null(borrow) & ay > iy) {
    # COMPUTE amount
    borrowing <- ctrl$value * borrow

    # STORE borrowing value in decision year
    track(tracking, "borrowing.isys", ay) <- borrowing

    # RECOVER borrowing
    borrowed <- tracking[[1]]["borrowing.isys", ac(ay - 1)]

    # CORRECT for borrowed
    ctrl$value <- ctrl$value - ifelse(is.na(borrowed), 0, borrowed)

    # APPLY borrowing
    ctrl$value <- ctrl$value + borrowing
  }

  # CALL split.is
  res <- split.is(stk, ctrl, split, quant = "catch", args, tracking)

  return(list(ctrl = res$ctrl, tracking = res$tracking))
}
# }}}

# bank_borrow.is {{{

bank_borrow.is <- function(
    stk, ctrl, args, split, rate = NULL, diff = 0.15,
    tracking) {
  # DIMS
  spread(args)

  # STOP if frq > 1
  if (frq > 1) {
    stop("Banking & borrowing currently assumes annual management (frq=1)")
  }

  # projection years
  pys <- ac(seq(ay + management_lag, ay + management_lag + frq - 1))

  # GET TAC
  tac <- ctrl$value

  # TODO: CHECK status

  # CORRECT for previous borrowing or banking
  if (ay == iy) {
    pre <- c(areaSums(unitSums(seasonSums(window(catch(stk),
      start = ay - management_lag, end = ay - management_lag
    )))))

    # START tracking
    track(tracking, "borrowing.isys", ac(ay)) <- 0
    track(tracking, "banking.isys", ac(ay)) <- 0
  } else {
    pre <- c(tracking[[1]]["hcr", ac(ay)])

    # GET banked or borrowed amounts
    borrowed <- c(tracking[[1]]["borrowing.isys", ac(ay)])
    banked <- c(tracking[[1]]["banking.isys", ac(ay)])

    tac <- tac - borrowed + banked
  }

  # IF lower TAC, THEN borrow
  borrowing <- ifelse(tac < pre * (1 - diff), tac * rate, 0)

  # TRACK amount being borrowed from pys
  track(tracking, "borrowing.isys", pys) <- borrowing

  # IF higher TAC, THEN bank
  banking <- ifelse(tac > pre * (1 + diff), tac * rate, 0)

  # TRACK amount being banked into pys
  track(tracking, "banking.isys", pys) <- banking

  # CORRECT TAC
  ctrl$value <- tac + borrowing - banking

  # split
  res <- split.is(stk, ctrl, split, quant = "catch", args, tracking)

  return(list(ctrl = res$ctrl, tracking = res$tracking))
}
# }}}

# catch.sel<- {{{
setReplaceMethod(
  "catch.sel", signature(object = "FLFishery", value = "FLQuant"),
  function(object, catch = 1, value) {
    if (length(object) > 1 & missing(catch)) {
      warnming("catch.sel modified for first FLCatch in object")
    }

    catch.sel(object[[catch]]) <- value

    return(object)
  }
)

setReplaceMethod(
  "catch.sel", signature(object = "FLFisheries", value = "FLQuants"),
  function(object, catch = 1, value) {
    object <- FLFisheries(Map(function(x, y, z) {
      catch.sel(x[[z]]) <- y
      return(x)
    }, x = object, y = value, z = catch))

    return(object)
  }
)
# }}}

# sampling.oem {{{

#' sampling.oem
#'
#' Samples from an operating model to obtain catch, biology and abundance data
#'
#' This observation error model (OEM) function mimics the most common data
#' collection regime, in which catch-at-age and biology is sampled from the
#' population, and one or more indices of abundance are derived from surveys
#' or CPUE data.
#'
#' The FLStock object passed to *sampling.oem* by the *mp* function is
#' simplified to match the dimensions of that present in the *observations*
#' slot.
#'
#' @param stk An FLStock object as obtained by the call to *stock(om)*.
#' @param deviances A named list of deviances, see Details.
#' @param observations A named list of observations, see Details.
#' @param args Options and arguments passed on by *mp()*.
#' @param tracking The tracking object.
#' @return A named list with elements *stk*, *idx*, *observations* and *tracking*.
#' @author Iago Mosqueira (WUR) & Ernesto Jardim (MSC).
#' @seealso \link{mp}
#' @keywords function
#' @examples
#' data(sol274)
#' # Generate samples from year 2000:2016
#' sampling.oem(stock(om),
#'   deviances = deviances(oem),
#'   observations = observations(oem),
#'   args = list(y0 = 2000, dy = 2021, frq = 1), tracking = FLQuant()
#' )
sampling.oem <- function(
    stk, deviances, observations, stability = 1,
    wts = TRUE, update = FALSE, args, tracking) {
  # DIMENSIONS
  y0 <- ac(args$y0)
  dy <- ac(args$dy)
  dyrs <- ac(seq(args$dy - args$frq + 1, args$dy))

  if (!update) {
    dyrs <- dyrs[an(dyrs) > args$iy]
  }

  # CHECK inputs
  if (any(!c("stk", "idx") %in% names(observations))) {
    stop("observations(oem) must have elements 'stk' and 'idx'.")
  }
  # if(!any(!c("stk", "idx") %in% names(deviances)))
  #  stop("deviances(oem) must have elements 'stk' or 'idx'.")

  # TODO: GENERATE length samples from OM: catch and indices, BUT needs selex
  # - invALK
  # - lenSamples

  # SUBSET year range
  stk <- window(stk, start = y0, end = dy, extend = FALSE)
  obs <- window(observations$stk, start = y0, end = dy, extend = FALSE)
  idx <- lapply(observations$idx, function(x) {
    if (dims(x)$maxyear >= dy) {
      window(x, end = dy)
    } else {
      x
    }
  })

  # RETURN IF no new years and update=FALSE
  if (length(dyrs) == 0) {
    return(list(stk = obs, idx = idx, observations = observations, tracking = tracking))
  }

  # --- STK

  if (!is.null(deviances$stk)) {
    # APPLY deviances and ASSIGN to stk slots in dyrs
    for (i in names(deviances$stk)) {
      slot(stk, i)[, dyrs] <-
        do.call(i, list(object = stk))[, dyrs] %*% deviances$stk[[i]][, dyrs] + 1e-8
    }

    # COMPUTE aggregated slots
    landings(stk)[, dyrs] <- computeLandings(stk[, dyrs])
    discards(stk)[, dyrs] <- computeDiscards(stk[, dyrs])
    catch(stk)[, dyrs] <- computeCatch(stk[, dyrs])
  }

  # --- IDX

  # CHOOSE indices to be updated (maxyear >= dy)
  upi <- unlist(lapply(idx, function(x) unname(dims(x)$maxyear) >= args$dy))

  if (is.null(deviances$idx) | length(deviances$idx) == 0) {
    deviances$idx <- lapply(observations$idx, function(x) index.q(x) %=% 1)
  }

  # APPLY survey() to stk + index (x), devs.q (y), stability(z)
  idx[upi] <- Map(
    function(x, y, z) {
      dyrs <- intersect(dyrs, dimnames(y)$year)

      # CREATE survey obs
      res <- survey(stk[, dyrs], x[, dyrs],
        sel = sel.pattern(x)[, dyrs],
        index.q = index.q(x)[, dyrs] %*% y[, dyrs], stability = z
      )

      # ENSURE no zeroes coming, maybe from high Fs
      if (sum(index(res)[, dyrs]) == 0) {
        index(res)[, dyrs] <- sqrt(.Machine$double.eps)
      }

      # SET 0s to min / 2
      index(res)[index(res) == 0] <- c(min(index(res)[index(res) > 0] / 2))

      # ASSIGN observations
      # TODO: ONLY if index not available
      x[, dyrs] <- res

      return(window(x, end = dy))
    },
    x = idx[upi], y = deviances$idx[upi],
    z = rep(stability, length = length(idx))[upi]
  )


  # UPDATE observations for dyrs
  for (i in seq(idx[upi])) {
    observations$idx[upi][[i]][, dyrs] <- idx[upi][[i]][, dyrs]
  }

  # CHECK dimensions to simplify, on catch.n for multi-fleet FLStock
  simp <- dim(catch.n(observations$stk))[c(3, 4, 5)] !=
    dim(catch.n(stk))[c(3, 4, 5)]

  # SIMPLIFY stk to match dimensions of observations$stk
  if (any(simp)) {
    stk <- simplify(stk, c("unit", "season", "area")[simp], harvest = FALSE)
  }

  # UPDATE observations
  slots <- c(
    "landings", "discards", "catch", "landings.n", "discards.n",
    "catch.n", "landings.wt", "discards.wt", "catch.wt", "stock.wt"
  )

  # UPDATE wts or only catches?
  if (!wts) slots <- slots[1:6]

  for (i in slots) {
    slot(obs, i)[, dyrs] <- slot(stk, i)[, dyrs]
  }

  # STORE in OEM observations
  observations$stk[, dyrs] <- obs[, dyrs]

  # RETURN stk from obs, only update observed 'slots'
  list(stk = obs, idx = idx, observations = observations, tracking = tracking)
} # }}}

# cpues.ind {{{

cpues.ind <- function(stk, idx, nyears = 5, ayears = 3, index = 1, args, tracking) {
  # CALL cpue.ind
  mets <- lapply(setNames(index, nm = names(idx[index])), function(x) {
    cpue.ind(stk, idx, nyears = nyears, ayears = ayears, index = x, args, tracking)
  })

  ind <- lapply(mets, function(x) x$ind$index)

  return(list(stk = stk, ind = ind, tracking = mets[[1]]$tracking))
} # }}}

# meta.hcr {{{

meta.hcr <- function(
    stk, ind, target, ..., args, tracking,
    combination = function(x, y) (x + y) / 2) {
  # DIMS
  dy <- ac(args$dy)

  # HCRs arguments
  rargs <- list(...)

  # CHECK rargs have 'method'
  # TODO: if(!all(unlist(lapply(rargs, function(x) is(method(x), "function")))))
  #  stop("args must contain list with an element called 'method'")

  # CHECK list has list w/ method + matching args
  # if(!all(unlist(lapply(rargs, function(x) all(names(x)[!grepl("method",
  #   names(x))] %in% names(formals(x$method)))))))
  #   stop("elements in each args list must match arguments in hcr function")

  # APPLY target to each hcr
  rargs <- lapply(rargs, function(x) {
    args(x)$target <- target
    return(x)
  })

  # APPLY each hcr, returns list(list('ctrl', 'vector'))
  decs <- lapply(rargs, function(x) {
    do.call(x@method, c(list(
      ind = ind, stk = stk,
      args = args, tracking = tracking
    ), x@args))
  })

  # if(!Reduce(all.equal, lapply(decs, function(x) as.character(x$ctrl$quant))))
  #  stop("Individual hcrs must output the same quant (e.g. 'catch')")

  # COMBINE ctrls
  ctrl <- decs[[1]]$ctrl
  ctrl$value <- Reduce("+", lapply(decs, function(x) x[[1]]$value)) / length(decs)

  # TRACK
  track(tracking, "tac.hcr", dy) <- ctrl$value

  # return
  list(ctrl = ctrl, tracking = tracking)
}
# }}}


# get_tuned_mps  {{{
get_tuned_mps <- function(dt, pattern = "tune", select_name = "mean(C)", buffer_filter = c("cpue3", "cpue36", "cpue367")) {
  dt[, year := as.integer(year)]
  dt[, iter := as.integer(iter)]
  dt[, data := as.numeric(data)]

  # Filter to relevant MPs and metric
  tuned_dt <- dt[grepl(pattern, mp) & name == select_name]

  # Extract tuning level
  tuned_dt[, tune_level := fifelse(
    grepl("tune05", mp), "tune05",
    fifelse(
      grepl("tune06", mp), "tune06",
      fifelse(grepl("tune07", mp), "tune07", NA_character_)
    )
  )]

  # Extract buffer type
  tuned_dt[, buffer_type := fifelse(
    grepl("cpue36_buffer", mp), "cpue36",
    fifelse(
      grepl("cpue3_buffer", mp), "cpue3",
      fifelse(grepl("cpue367_buffer", mp), "cpue367", NA_character_)
    )
  )]

  tuned_dt[, tune_level := factor(tune_level, levels = c("tune05", "tune06", "tune07"))]
  tuned_dt[, buffer_type := factor(buffer_type, levels = c("cpue3", "cpue36", "cpue367"))]

  # Apply buffer_filter
  tuned_dt <- tuned_dt[buffer_type %in% buffer_filter]

  return(tuned_dt)
}
# }}}


# cpuescore.ind {{{
# Calculated over mean and sd on refyrs

cpuescore.ind <- function(stk, idx, index = 1, refyrs = NULL, args, tracking) {
  # ARGS
  ay <- args$ay
  dlag <- args$data_lag
  dy <- ay - dlag
  # TODO: CHECK for frq>1

  # GET metric until dy
  met <- window(biomass(idx[[index]])[1, ], end = dy)

  if (!is.null(refyrs)) {
    ref <- met[, ac(refyrs)]
  } else {
    ref <- met
  }

  ind <- FLQuants(zscore = (met[, ac(dy)] %-% yearMeans(ref)) %/%
    sqrt(yearVars(ref)))

  track(tracking, "cpue.ind", ac(ay)) <- ind$zscore

  return(list(stk = stk, ind = ind, tracking = tracking, cpue = met))
}
# }}}

# cpuescore2.ind {{{
# Calculated over mean and sd on passed as arguments

cpuescore2.ind <- function(stk, idx, index = 1, mean, sd, args, tracking) {
  # ARGS
  ay <- args$ay
  dlag <- args$data_lag
  dy <- ay - dlag
  # TODO: CHECK for frq>1

  # GET metric until dy
  met <- biomass(idx[[index]])[, ac(dy)]

  ind <- FLQuants(zscore = (met[, ac(dy)] %-% mean) %/% sd)

  track(tracking, "cpue.ind", ac(ay)) <- ind$zscore

  return(list(stk = stk, ind = ind, tracking = tracking, cpue = met))
}
# }}}

# cpuescores.ind {{{

cpuescores.ind <- function(
    stk, idx, index = 1, dlag = rep(args$data_lag, length(index)),
    refyrs = NULL, args, tracking) {
  # ARGS
  dlag <- setNames(dlag, nm = names(idx)[index])
  ay <- args$ay
  dy <- ay - dlag
  # TODO: CHECK for frq>1

  res <- as.list(setNames(nm = names(idx)[index]))

  # GET metric until dy
  for (i in names(res)) {
    met <- window(biomass(idx[[i]])[1, ], end = dy[i])

    if (!is.null(refyrs)) {
      ref <- met[, ac(refyrs)]
    } else {
      ref <- met
    }

    res[[i]] <- (met[, ac(dy[i])] %-% yearMeans(ref)) %/% sqrt(yearVars(ref))

    dimnames(res[[i]])$year <- max(dy)

    track(tracking, paste0("score.mean.", i), ac(ay)) <- yearMeans(ref)
    track(tracking, paste0("score.sd.", i), ac(ay)) <- sqrt(yearVars(ref))
    track(tracking, paste0("score.ind.", i), ac(ay)) <- res[[i]]
  }

  return(list(stk = stk, ind = res, tracking = tracking, cpue = met))
}

# }}}

zscore <- function(x) (x %-% yearMeans(x)) %/% sqrt(yearVars(x))

cpuescorejim.ind <- function(stk, idx, index = 1, mean, sd, args, tracking) {
  # ARGS
  ay <- args$ay
  dlag <- args$data_lag
  dy <- ay - dlag
  # TODO: CHECK for frq>1

  # GET metric until dy
  met <- biomass(idx[[index]])[, ac(dy)]

  ind <- FLQuants(zscore = (met[, ac(dy)] %-% mean) %/% sd)
  ind$zscore <- ifelse(  ind$zscore > 3, 3, 
                         ifelse(  ind$zscore < -3, -3, ind$zscore))

  track(tracking, "cpue.ind", ac(ay)) <- ind$zscore

  return(list(stk = stk, ind = ind, tracking = tracking, cpue = met))
}
# }}}


# periodsPerformance {{{

periodsPerformance <- function(x, periods) {
  # COERCE tio list
  periods <- as.list(periods)

  years <- unlist(lapply(periods, function(x) {
    if (length(x) > 1) {
      paste(x[1], substr(rev(x)[1], 3, 4), sep = "-")
    } else {
      x
    }
  }))

  # ASSIGN names if missing
  if (is.null(names(periods))) {
    names(periods) <- years
  }

  res <- rbindlist(Map(
    function(pe, na, ye) {
      x[year %in% pe, .(data = mean(data, na.rm = TRUE), period = na, year = ye),
        by = .(type, mp, statistic, name, desc, iter)
      ]
    },
    pe = periods, na = names(periods), ye = years
  ))

  return(res)
}
# }}}

# labelPerformance {{{

labelPerformance <- function(dat, labels = NULL) {
  # NO label, use mp
  if (is.null(labels)) {
    dat[, label := mp]
    return(dat[])
    # 'numeric', set as sequence in sort order
  } else if (identical(labels, "numeric")) {
    labels <- data.table(mp = sort(unique(dat$mp)), label = seq(unique(dat$mp)))
    # VECTOR, assign names by sort(unique)
  } else if (is.vector(labels)) {
    labels <- data.table(mp = sort(unique(dat$mp)), label = labels)
    # SET as data.table JIC
  } else {
    labels <- data.table(labels)
  }

  # CHECK dims
  if (!all(unique(dat[, mp]) %in% unique(labels[, mp]))) {
    stop("'mp' names in both tables do not match")
  }

  # MERGE by mp
  dat <- merge(dat[, !"label"], labels[, .(mp, label)], by = "mp")

  return(dat[])
}

