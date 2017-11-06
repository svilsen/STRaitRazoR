if (FALSE) {
    require("readr", "dplyr", "tidyr")

    input <- "~/MPS/R705-A501_S1_L001_R1_001.fastq"
    outputLocation <- "~/allSeqs.txt"

    STRaitRazorLinux(input, output, configFile = "forenseq")
    tt <- STRaitRazorLinuxSTRMPS(output)
}
