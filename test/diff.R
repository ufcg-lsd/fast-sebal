suscess <- require(raster)

if(suscess) {

    options(digits = 10)

    args = commandArgs(trailingOnly=TRUE)

    path.tiff.r <- args[1]
    path.tiff.c <- args[2]

    if (length(args) == 3) {
       sink(args[3])
    } else {
       sink("out")
    }
    
    TiffR <- raster(path.tiff.r)
    TiffC <- raster(path.tiff.c)

    cat("----------- R interval -----------", "\n")
    maxR <- max(TiffR, na.rm = TRUE)
    minR <- min(TiffR, na.rm = TRUE)
    toPrint <- paste("[", maxR, ",", minR, "]")
    cat(toPrint, "\n")

    cat("----------- C++ interval -----------", "\n")
    maxC <- max(TiffC, na.rm = TRUE)
    minC <- min(TiffC, na.rm = TRUE)
    toPrint <- paste("[", maxC, ",", minC, "]")
    cat(toPrint, "\n")

    cat("----------- Absolute error -----------", "\n")
    TiffDiff <- abs(TiffR[] - TiffC[])
    maxDiff <- max(TiffDiff, na.rm = TRUE)
    cat(paste("Max absolute error:", maxDiff), "\n")

    cat("----------- Percentual error -----------", "\n")
    TiffPer <- (TiffDiff / abs(TiffR[])) * 100
    maxDiff <- max(TiffPer, na.rm = TRUE)

    if (maxDiff == Inf) {
       
        TiffR[TiffR == 0] <- NA
        TiffPer <- (TiffDiff / abs(TiffR[])) * 100
        maxDiff <- max(TiffPer, na.rm = TRUE)

    }

    cat(paste("Max percentual error:", maxDiff), "\n")

    distribution <- c(0, 0, 0, 0, 0)
    names(distribution) <- c("equals 0", "0% ~ 0.1%", "0.1% ~ 1%", "1% ~ 10%", "bigger than 10%")
    distribution[1] <- length(which(TiffPer[] == 0))
    distribution[2] <- length(which(TiffPer[] > 0 && TiffPer[] <= 0.1))
    distribution[3] <- length(which(TiffPer[] > 0.1 && TiffPer[] <= 1))
    distribution[4] <- length(which(TiffPer[] > 1 && TiffPer[] <= 10))
    distribution[5] <- length(which(TiffPer[] > 10))

    cat(distribution, "\n")

    cat("----------- Percentual error ignoring values less than 1e-4 -----------", "\n")

    TiffC[TiffC < 1e-4] <- NA
    TiffR[TiffR < 1e-4] <- NA
    TiffDiff <- abs(TiffR[] - TiffC[])
    TiffPer <- (TiffDiff / abs(TiffR[])) * 100
    maxDiff <- max(TiffPer, na.rm = TRUE)

    cat(paste("Max percentual error:", maxDiff), "\n")

    distribution <- c(0, 0, 0, 0, 0)
    names(distribution) <- c("equals 0", "0% ~ 0.1%", "0.1% ~ 1%", "1% ~ 10%", "bigger than 10%")
    distribution[1] <- length(which(TiffPer[] == 0))
    distribution[2] <- length(which(TiffPer[] > 0 && TiffPer[] <= 0.1))
    distribution[3] <- length(which(TiffPer[] > 0.1 && TiffPer[] <= 1))
    distribution[4] <- length(which(TiffPer[] > 1 && TiffPer[] <= 10))
    distribution[5] <- length(which(TiffPer[] > 10))

    cat(distribution, "\n")

} else {

    print("raster package not installed")

}