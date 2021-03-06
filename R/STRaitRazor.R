utils::globalVariables(c("Allele", "C1", "C2", "Coverage", "Length", "Marker", "MarkerAllele", "MotifLength",
                         "Region", "Type"))

#' @title STRaitRazor v3
#'
#' @description Simple wrap function for calling the command line tool STRaitRazor v3 directly from R (on linux systems).
#'
#' @param inputLocation Path to fastq file.
#' @param outputLocation Path to the .txt returned by the STRait Razor v3 tool
#' @param configFile The STRait Razor '.config' file to be used in the analysis (only ForenSeq and PowerSeq available).
#' @param commandLineArguments Any command line arguments passed to STRait Razor excluding -c and -p.
#' @param numberOfThreads The number of threads used to run the analysis.
STRaitRazorLinux <- function(inputLocation, outputLocation, commandLineArguments = NULL, configFile  = NULL, numberOfThreads = 4) {
    configFile = ifelse(is.null(configFile), "forenseq", configFile)
    config = switch(tolower(configFile),
                    "forenseq" = "Forenseq.config",
                    "powerseq" = "Powerseq.config")

    STRaitPath <- system.file("bin/STRaitRazor", "", package = "STRaitRazoR")
    if (length(STRaitPath) == 0) {
        stop("Couldn't find STRait Razor commandline tool.")
    }

    commandLineArguments = ifelse(is.null(commandLineArguments), "-c", commandLineArguments)
    if (numberOfThreads > 1) {
        commandLineArguments <- paste("-p", numberOfThreads, commandLineArguments)
    }

    if (tolower(Sys.info()['sysname']) == "linux") {
        commandLineSTRaitRazor <- paste(STRaitPath, "str8rzr ", commandLineArguments, " ", STRaitPath, config, " ", inputLocation, " > ", outputLocation, sep = "")
        isWindows <- FALSE
    }
    else {
        stop("Not implemented for OSX or windows.")
    }

    system(commandLineSTRaitRazor)
}

#' @title STRaitRazor to STRMPS control function
#'
#' @description Contains the arguments passed to the \link{STRaitRazorLinux} function, when it is called through the \link{STRaitRazorSTRMPS} function.
#'
#' @param commandLineArguments Command line arguments passed to STRaitRazor v3.
#' @param configFile The STRaitRazor config file used for flanking regoin identification.
#' @param numberOfThreads The total number of cores used for parallelisation.
#'
#' @return A list containing default input.
STRaitRazorSTRMPS.control <- function(commandLineArguments = NULL, configFile  = NULL, numberOfThreads = 4) {
    return(list(commandLineArguments = commandLineArguments, configFile = configFile, numberOfThreads = numberOfThreads))
}

#' A string coverage list
#'
#' A list of tibbles, one for every marker, used to contain the sequencing information of STR MPS data.
setClass("stringCoverageList")

#' @title STRaitRazor v3 output to STRMPS data structure.
#' @description A function for converting the STRaitRazor output to an object useable by the STRMPS-package.
#'
#' @param inputLocation Path to fastq file.
#' @param outputLocation Path to stored output. If NULL a temp folder is used and the \link{STRaitRazorLinux} function is called.
#' @param control A control object containing additional input for the \link{STRaitRazorLinux} function.
#'
#' @return A list of tibbles used by the STRMPS package.
STRaitRazorSTRMPS <- function(inputLocation, outputLocation = NULL, control = STRaitRazorSTRMPS.control()) {
    STRaitPath <- system.file("bin/STRaitRazor", "", package = "STRaitRazoR")
    if (STRaitPath == "") {
        stop("Couldn't find STRaitRazor commandline tool.")
    }

    if (is.null(outputLocation)) {
        outputLocation <- paste(tempdir(), "all_seqs_temp.txt", sep = "/")
        STRaitRazorLinux(inputLocation, outputLocation = outputLocation, commandLineArguments = control$commandLineArguments,
                    configFile = control$configFile, numberOfThreads = control$numberOfThreads)
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


