#!/bin/bash

# Default vars
ECHO=$(echo [$0])
TYPE="qdnaseq"
NOWGD=""
TCN=""
CORES=8
DEBUG=FALSE
ARGS=""

# Stop/error function
stopFunction(){
    helpFunction
    exit 1
}

# Termination trap
trap 'stopFunction' SIGINT SIGTERM

# Help function
helpFunction(){
    echo -e "${ECHO}[`date "+%H:%M:%S"`][INIT] Wrapper script for multi-sample and variable input types for medicc2"
    echo -e "Options:"
    echo -e "\tArguments:"
    echo -e "\t\t-t\tinput type (either qdnaseq | segment | medicc)"
    echo -e "\t\t-i\tFull path to input file matching input type specified by -t"
    echo -e "\t\t-o\tFull path output directory for Medicc2 results"
    echo -e "\t\t-m\tFull path to meta file detailing patient-sample relationship, where -t is not medicc"
    echo -e "\t\t-a\tAdditional arguments passed to medicc2 (must be quoted)"
    echo -e "\tFlags:"
    echo -e "\t\t-w\tSet medicc2 to not use whole genome duplication in tree calculations"
    echo -e "\t\t-n\tSet midicc2 to use total copy number (default column name 'cn_a')"
    echo -e "\t\t-wn\tSet both options above"
    echo -e ""
}

makeDirs(){
    if ! [ -d "${OUTPUT}input_files/" ]; then
        mkdir ${OUTPUT}input_files/
    fi

    if ! [ -d "${OUTPUT}medicc2_output/" ]; then
        mkdir ${OUTPUT}medicc2_output/
    fi
}

set -e

if [[ $# -eq 0 ]]; then
	helpFunction; exit 0
fi

while getopts "hwdni:o:t:c:m:a:" opt; do
  case $opt in
    h ) helpFunction ; exit 0
    ;;
    i ) INPUT="$OPTARG"
    ;;
    o ) OUTPUT="$OPTARG"
    ;;
    t ) TYPE="$OPTARG"
	[[ ! $TYPE =~ qdnaseq|segment|medicc ]] && {
            echo -e "${ECHO}[`date "+%H:%M:%S"`][ERROR] Incorrect input type provided"
            stopFunction
        }
    ;;
    m ) META="$OPTARG"
    ;;
    c ) CORES="$OPTARG"
    ;;
    w ) NOWGD=" --no-wgd"
    ;;
    n ) TCN=" --total-copy-number -a cn_a"
    ;;
    a ) ARGS="$OPTARG"
    ;;
    d ) DEBUG=TRUE
    ;;
    :  ) echo -e "${ECHO}[`date "+%H:%M:%S"`][ERROR] Invalid option: $OPTARG requires an argument" ; stopFunction
    ;;
    #\? ) echo "Invalid option: $OPTARG" ; stopFunction
    #;;
    #*  ) echo -e "${ECHO}[`date "+%H:%M:%S"`][ERROR] Unimplemented option: $OPTARG" ; stopFunction
    #;;
  esac
done

## VAR PRINT CHECK
if [ "$DEBUG" == TRUE ]; then
    echo -e "\${INPUT}=${INPUT}"
    echo -e "\${OUPUT}=${OUTPUT}"
    echo -e "\${META}=${META}"
    echo -e "\${TYPE}=${TYPE}"
    echo -e "\${NOWGD}=${NOWGD}"
    echo -e "\${TCN}=${TCN}"
    echo -e "\${DEBUG}=${DEBUG}"
    echo -e "\${CORES}=${CORES}"
    echo -e "\${ARGS}=${ARGS}"
    exit 0
fi

## Check input file
if ! [ -f "${INPUT}" ]; then
    echo -e "${ECHO}[`date "+%H:%M:%S"`][ERROR] Input file does not exist"
    stopFunction
elif ! [ -r "${INPUT}" ]; then
    echo -e "${ECHO}[`date "+%H:%M:%S"`][ERROR] Input file not readable"
    stopFunction
fi

## Check META file
if [ "$TYPE" != "medicc" ]; then
    if ! [ -f "${META}" ]; then
        echo -e "${ECHO}[`date "+%H:%M:%S"`][ERROR] Metadata file does not exist"
        stopFunction
    elif ! [ -r "${META}" ]; then
        echo -e "${ECHO}[`date "+%H:%M:%S"`][ERROR] Metadata file not readable"
        stopFunction
    fi
fi

## Check output directory exists and writable 
if ! [ -d "${OUTPUT}" ]; then
    echo -e "${ECHO}[`date "+%H:%M:%S"`][ERROR] Output directory does not exist"
    stopFunction
elif ! [ -w "${OUTPUT}" ]; then
    echo -e "${ECHO}[`date "+%H:%M:%S"`][ERROR] Output directory not writable"
    stopFunction
fi

## Check medicc2 is available
if ! [ -x "$(command -v medicc2)" ]; then
    echo -e "${ECHO}[`date "+%H:%M:%S"`][ERROR] Medicc2 unavailable - check conda environment"
    stopFunction
fi

echo -e "${ECHO}[`date "+%H:%M:%S"`][INIT] Making working directories"
makeDirs

if [ "$TYPE" == "qdnaseq" ]; then
    ## Check programmes available
    echo -e "${ECHO}[`date "+%H:%M:%S"`][INIT] Input type set to QDNASEQ object"
    if ! [ -x "$(command -v Rscript)" ]; then
        echo -e "${ECHO}[`date "+%H:%M:%S"`][ERROR] Rscript unavailable - check conda environment"
        stopFunction
    fi
    echo -e "${ECHO}[`date "+%H:%M:%S"`][MAIN] Extracting medicc2 input format"
    ## Build output directory structure
    Rscript scripts/qdnaseq_to_medicc_format_script.R ${INPUT} ${OUTPUT}input_files/ ${META}
elif [ "$TYPE" == "segment" ]; then
    echo -e "${ECHO}[`date "+%H:%M:%S"`][INIT] Input type set to segment table"
    if ! [ -x "$(command -v Rscript)" ]; then
        echo -e "${ECHO}[`date "+%H:%M:%S"`][ERROR] Rscript unavailable - check conda environment"
        stopFunction
    fi
    echo -e "${ECHO}[`date "+%H:%M:%S"`][MAIN] Extracting medicc2 input format"
    ## Build output directory structure
    Rscript scripts/segtable_to_medicc_format_script.R ${INPUT} ${OUTPUT}input_files/ ${META}
elif [ "$TYPE" == "medicc" ]; then
    echo -e "${ECHO}[`date "+%H:%M:%S"`][INIT] Input type set to medicc"
    cp ${INPUT} ${OUTPUT}input_files/
fi

echo -e "${ECHO}[`date "+%H:%M:%S"`][MAIN] Running medicc2"
find ${OUTPUT}input_files/ -name "*" -type f | xargs -P ${CORES}  -n 1 -I {} medicc2 {} ${OUTPUT}medicc2_output/${NOWGD}${TCN} ${ARGS}

echo -e "${ECHO}[`date "+%H:%M:%S"`][MAIN] Generating summary information"
## Combine trees
cat ${OUTPUT}medicc2_output/*final_tree.new > ${OUTPUT}medicc2_output/all.trees.new

echo -e "${ECHO}[`date "+%H:%M:%S"`][MAIN] Pipeline complete"
## END
