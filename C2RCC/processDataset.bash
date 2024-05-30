#!/bin/bash
# enable next line for debugging purpose
# set -x 

# got this script from this website:
# https://senbox.atlassian.net/wiki/spaces/SNAP/pages/70503475/Bulk+Processing+with+GPT

############################################
# User Configuration
############################################

# adapt this path to your needs
export PATH=/Users/daa5@stir.ac.uk/Applications/esa-snap/bin:$PATH
gptPath="gpt"

############################################
# Command line handling
# these parameters are how you call this script from the terminal
# eg. processDataset.bash resample_s2.xml resample_20m.properties "/Eodata/toProcess" "/Eodata/toProcess/output" resampled20m
############################################

# first parameter is a path to the graph xml
graphXmlPath="$1"

# second parameter is a path to a parameter file
parameterFilePath="$2"

# use third parameter for path to source products
sourceDirectory="$3"

# use fourth parameter for path to target products
targetDirectory="$4"

# the fifth parameter is a file prefix for the target product name, typically indicating the type of processing
targetFilePrefix="$5"

   
############################################
# Helper functions
############################################
removeExtension() {
    file="$1"
    echo "$(echo "$file" | sed -r 's/\.[^\.]*$//')"
}


############################################
# Main processing
############################################

# Create the target directory
mkdir -p "${targetDirectory}"

# the d option limits the elemeents to loop over to directories. Remove it, if you want to use files.
for F in $(ls -1d "${sourceDirectory}"/S2*.SAFE); do
  sourceFile="$(realpath "$F")"
  targetFile="${targetDirectory}/${targetFilePrefix}_$(removeExtension "$(basename ${F})").nc"
  ${gptPath} ${graphXmlPath} -e -p ${parameterFilePath} -t ${targetFile} ${sourceFile}
done
