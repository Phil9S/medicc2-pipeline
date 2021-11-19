#!/bin/bash

# Default vars
ECHO=$(echo [$0])
TYPE="qdnaseq"

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
}

makeDirs(){
    if ! [ -d "${OUPUT}input_files/" ]; then
        mkdir ${OUTPUT}input_files/
    fi

    if ! [ -d "${OUPUT}medicc2_output/" ]; then
        mkdir ${OUTPUT}medicc2_output/
    fi
}


set -e

if [[ $# -eq 0 ]]; then
	helpFunction; exit 0
fi

while getopts "hwci:o:t:" opt; do
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
    : ) echo -e "${ECHO}[`date "+%H:%M:%S"`][ERROR] Invalid option: $OPTARG requires an argument" ; stopFunction
    ;;
    #\? ) echo "Invalid option: $OPTARG" ; stopFunction
    #;;
    #*  ) echo "Unimplemented option: $OPTARG" ; stopFunction
    #;;
  esac
done

## Check input file
if ! [ -f "${INPUT}" ]; then
    echo -e "${ECHO}[`date "+%H:%M:%S"`][ERROR] Input file does not exist"
    stopFunction
elif ! [ -r "${INPUT}" ]; then
    echo -e "${ECHO}[`date "+%H:%M:%S"`][ERROR] Input file not readable"
    stopFunction
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

if [ "$TYPE" == "qdnaseq" ]; then
    ## Check programmes available
    if ! [ -x "$(command -v Rscript)" ]; then
        echo -e "${ECHO}[`date "+%H:%M:%S"`][ERROR] Rscript unavailable - check conda environment"
        stopFunction
    fi
    echo -e "QDNASEQ"
    ## Build output directory structure
    makeDirs
    #Rscript scripts/qdnaseq_to_medicc_format_script.R ${INPUT} ${OUPUT}input_files/
fi

echo "END"
