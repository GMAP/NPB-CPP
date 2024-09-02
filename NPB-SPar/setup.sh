#!/bin/bash

# Defining the library directory
LIB_DIR="./lib"
SPAR_DIR="$LIB_DIR/spar"

# Checking if the library directory exists
if [ -d "$SPAR_DIR" ]; then
    echo "The spar folder already exists. No need to download it again."
else
    # Creating the lib directory if it doesn't exist
    if [ ! -d "$LIB_DIR" ]; then
        mkdir -p "$LIB_DIR"
    fi

    # Navigating to the lib directory
    cd "$LIB_DIR"

    # Cloning the SPar repository
    git clone https://github.com/GMAP/SPar

    # Renaming the cloned folder to spar
    mv SPar spar

    # Navigating to the spar directory
    cd spar

    # Compiling the library
    make
fi
