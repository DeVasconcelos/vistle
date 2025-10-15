#!/bin/bash

RUN_SCRIPT_DIR=$PWD

# check if at least two arguments are given
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <.vsl file> <prefix for output folder name>"
    exit 1
fi

VSL_MAP="$1"
OUTPUT_PREFIX="$2"

OUTPUT_DIRECTORY="${RUN_SCRIPT_DIR}/${OUTPUT_PREFIX}-n$GRID_SIZE-it$ITERATIONS"

# warn user that output directory will be overwritten if it exists
if [ -d "$OUTPUT_DIRECTORY" ]; then
    echo "Warning: Folder '$OUTPUT_DIRECTORY' already exists. Its contents will be overwritten by this script."
    read -p "Do you want to continue? (y/n): " answer
    case $answer in
        [Yy]*)
            echo "Continuing execution..."
            rm -rf $OUTPUT_DIRECTORY
            mkdir -p $OUTPUT_DIRECTORY;;
        [Nn]*)
            echo "Exiting script!"
            exit 1;;
        *)
            echo "Invalid answer. Exiting script!"
            exit 1;;
    esac
fi

cd $VISTLE_BUILD_TO_USE
# load environment so vistle is set up correctly (note that direnv is intended for interactive shells)
direnv allow . && eval "$(direnv export bash)" 

for i in $(seq "$ITERATIONS"); do
    export VISTLE_LOGFILE=$OUTPUT_DIRECTORY/log-iteration-$i.txt
    ./bin/vistle -b $VSL_MAP
done

# calculate average run times from log files
cd $RUN_SCRIPT_DIR
eval "$(direnv export bash)" # unload environment
python extract_times.py $OUTPUT_DIRECTORY

if [[ -f "times.txt" ]]; then
    echo "Now storing result into $OUTPUT_DIRECTORY..."
    mv times.txt $OUTPUT_DIRECTORY
else 
    echo "WARNING: extract_times.py did not produce output."
fi
