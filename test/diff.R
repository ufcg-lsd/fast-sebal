suscess <- require(raster)

if(suscess) {

    options(digits = 10)
    options(width = 200)

    args = commandArgs(trailingOnly=TRUE)

    path.tiff.r <- args[1]
    path.tiff.c <- args[2]

    if (length(args) == 3) {
       sink(args[3])
    } else {
       sink("out")
    }

    cat("----------- Tiff Informations -----------", "\n")
    
    TiffR <- raster(path.tiff.r)
    TiffC <- raster(path.tiff.c)
    number.cell <- ncell(TiffR)

    TiffC[is.nan(TiffR[])] <- NaN
    TiffR[is.nan(TiffC[])] <- NaN

    cat(paste("Number of cell:", number.cell), "\n")

    cat("\n")

    cat("----------- R interval -----------", "\n")
    maxR <- max(TiffR[], na.rm = TRUE)
    minR <- min(TiffR[], na.rm = TRUE)
    toPrint <- paste("[", minR, ",", maxR, "]")
    cat(toPrint, "\n")

    cat("\n")

    cat("----------- C++ interval -----------", "\n")
    maxC <- max(TiffC[], na.rm = TRUE)
    minC <- min(TiffC[], na.rm = TRUE)
    toPrint <- paste("[", minC, ",", maxC, "]")
    cat(toPrint, "\n")

    cat("\n")

    cat("----------- Absolute error -----------", "\n")
    TiffDiff <- abs(TiffR[] - TiffC[])
    maxDiff <- max(TiffDiff[], na.rm = TRUE)
    cat(paste("Max absolute error:", maxDiff), "\n")

    cat("\n")

    cat("----------- Percentual error -----------", "\n")
    TiffPer <- (TiffDiff[] / abs(TiffR[])) * 100
    maxDiff <- max(TiffPer[], na.rm = TRUE)

    if (maxDiff == Inf) {
       
        TiffR[TiffR[] == 0] <- NA
        TiffPer <- (TiffDiff[] / abs(TiffR[])) * 100
        maxDiff <- max(TiffPer[], na.rm = TRUE)

    }

    cat(paste("Max percentual error:", maxDiff), "\n")

    distribution <- c(0, 0, 0, 0, 0, 0)
    names(distribution) <- c("equals 0", "0% ~ 0.1%", "0.1% ~ 1%", "1% ~ 10%", "bigger than 10%", "NA")
    distribution[1] <- length(TiffPer[!is.na(TiffPer[]) & TiffPer[] == 0])
    distribution[2] <- length(TiffPer[!is.na(TiffPer[]) & TiffPer[] > 0 & TiffPer[] <= 0.1])
    distribution[3] <- length(TiffPer[!is.na(TiffPer[]) & TiffPer[] > 0.1 & TiffPer[] <= 1])
    distribution[4] <- length(TiffPer[!is.na(TiffPer[]) & TiffPer[] > 1 & TiffPer[] <= 10])
    distribution[5] <- length(TiffPer[!is.na(TiffPer[]) & TiffPer[] > 10])
    distribution[6] <- length(TiffPer[is.na(TiffPer[])])

    if (sum(distribution) != number.cell) {
       print("Something wrong here!")
    }

    distribution <- distribution * 100 / number.cell

    options(digits = 5)

    print(distribution)

    options(digits = 10)

    cat("\n")

    cat("----------- Percentual error ignoring values less than 1e-4 -----------", "\n")

    TiffC[TiffC[] < 1e-4] <- NA
    TiffR[TiffR[] < 1e-4] <- NA
    TiffDiff <- abs(TiffR[] - TiffC[])
    TiffPer <- (TiffDiff[] / abs(TiffR[])) * 100
    maxDiff <- max(TiffPer[], na.rm = TRUE)

    cat(paste("Max percentual error:", maxDiff), "\n")

    distribution <- c(0, 0, 0, 0, 0, 0)
    names(distribution) <- c("equals 0", "0% ~ 0.1%", "0.1% ~ 1%", "1% ~ 10%", "bigger than 10%", "NA")
    distribution[1] <- length(TiffPer[!is.na(TiffPer[]) & TiffPer[] == 0])
    distribution[2] <- length(TiffPer[!is.na(TiffPer[]) & TiffPer[] > 0 & TiffPer[] <= 0.1])
    distribution[3] <- length(TiffPer[!is.na(TiffPer[]) & TiffPer[] > 0.1 & TiffPer[] <= 1])
    distribution[4] <- length(TiffPer[!is.na(TiffPer[]) & TiffPer[] > 1 & TiffPer[] <= 10])
    distribution[5] <- length(TiffPer[!is.na(TiffPer[]) & TiffPer[] > 10])
    distribution[6] <- length(TiffPer[is.na(TiffPer[])])

    if (sum(distribution) != number.cell) {
       print("Something wrong here!")
    }

    distribution <- distribution * 100 / number.cell
    
    options(digits = 5)

    print(distribution)

} else {

    print("raster package not installed")

}