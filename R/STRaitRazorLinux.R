#' @title STRait Razor Linux
#'
#' @description Simple wrap function for calling the command line tool STRaitRazor directly from R in a linux OS.
#'
#' @param inputLocation Path to fastq file.
#' @param outputLocation Path to the .txt returned by the STRait Razor v3 tool
#' @param configFile The STRait Razor '.config' file to be used in the analysis (only ForenSeq and PowerSeq available).
#' @param commandLineArguments Any command line arguments passed to STRait Razor excluding -c and -p.
#' @param numberOfThreads The number of threads used to run the analysis.
#'
#' @param return Nothing is returned to R, but a .txt file is created and stored at the output location.
STRaitRazorLinux <- function(inputLocation, outputLocation, commandLineArguments = NULL, configFile  = NULL, numberOfThreads = 4) {
    configFile = ifelse(is.null(configFile), "forenseq", configFile)
    config = switch(tolower(configFile),
                    "forenseq" = "Forenseq.config",
                    "forenseqstrs" = "ForenseqSTRs.config",
                    "forenseqautosomal" = "ForenseqAUTOSOMAL.config",
                    "powerseq" = "Powerseq.config")
    libPaths <- .libPaths()

    STRaitPath <- paste(libPaths, "STRaitRazoR/STRaitRazor", sep = "/")[which(dir.exists(paste(libPaths, "STRaitRazoR/STRaitRazor", sep = "/")))]
    if (length(STRaitPath) == 0) {
        stop("Couldn't find STRait Razor commandline tool.")
    }

    commandLineArguments = ifelse(is.null(commandLineArguments), "-c", commandLineArguments)
    if (numberOfThreads > 1) {
        commandLineArguments <- paste("-p", numberOfThreads, commandLineArguments)
    }

    commandLineSTRaitRazor <- paste(STRaitPath, "/str8rzr ", commandLineArguments, " ", STRaitPath, "/", config, " ", inputLocation, " > ", outputLocation, sep = "")
    system(commandLineSTRaitRazor)
}


setClass("stringCoverageList")

#' @title STRait Razor output to STRMPS
#' @description A function for converting the output created by STRait Razor to an object useable by the STRMPS-package.
#'
#' @param outputLocation Location of the output created by the \link{STRaitRazorLinux}-function (or any run of the STRait Razor v3 command line tool).
#'
#' @return A list of tibbles used in the STRMPS package.
STRaitRazorLinuxSTRMPS <- function(outputLocation) {
    libPaths <- .libPaths()

    STRaitPath <- paste(libPaths, "STRaitRazoR/STRaitRazor", sep = "/")[which(dir.exists(paste(libPaths, "STRaitRazoR/STRaitRazor", sep = "/")))]
    if (length(STRaitPath) == 0) {
        stop("Couldn't find STRait Razor commandline tool.")
    }

    addedInformation <- bind_rows(read_delim(paste(STRaitPath, "/Forenseq.config", sep = ""), "\t", skip = 1,
                                   col_names = c("Marker", "Type", "F", "R", "M", "MotifLength", "O"),
                                   col_types = list(col_character(), col_character(), col_character(), col_character(), col_character(), col_integer(), col_integer())) %>%
        select(Marker, Type, MotifLength), read_delim(paste(STRaitPath, "/Powerseq.config", sep = ""), "\t", skip = 1,
               col_names = c("Marker", "Type", "F", "R", "M", "MotifLength", "O"),
               col_types = list(col_character(), col_character(), col_character(), col_character(), col_character(), col_integer(), col_integer())) %>%
        select(Marker, Type, MotifLength)) %>% distinct(Marker, .keep_all = TRUE)

    stringCoverageTibble <- read_delim(outputLocation, "\t",
                                       col_names = c("MarkerAllele", "Length", "Region", "C1", "C2"),
                                       col_types = list(col_character(), col_character(), col_character(), col_integer(), col_integer())) %>%
        mutate(Coverage = C1 + C2) %>% select(-C1, -C2) %>%
        separate(MarkerAllele, c("Marker", "Allele"), sep = ":") %>% separate(Length, c("Length", "D"), sep = " ") %>%
        left_join(addedInformation, by = "Marker") %>% mutate(Marker = toupper(Marker), Allele = as.numeric(Allele)) %>%
        select(Marker, Type, Allele, Region, MotifLength, Coverage) %>%
        arrange(Marker, Allele, Coverage)

    stringCoverageList <- split(stringCoverageTibble, as.factor(stringCoverageTibble$Marker))
    stringCoverageList <- stringCoverageList[order(sapply(stringCoverageList, function(x) unique(x$Type)))]

    class(stringCoverageList) <- "stringCoverageList"
    return(stringCoverageList)
}


