#!/bin/bash

# Check and create input and input/FASTQ folders
if [ ! -d "input" ] || [ ! -d "input/FASTQ" ]; then
    mkdir -p "input/FASTQ"
    echo "Created folder: input/FASTQ"
    mkdir -p "input"
    echo "Created folder: input"
fi

# Check and create output folder
if [ ! -d "output" ]; then
    mkdir -p "output"
    echo "Created folder: output"
fi

# Check and create resources folder
if [ ! -d "resources" ]; then
    mkdir -p "resources"
    echo "Created folder: resources"
fi

# Check and create scripts.dev folder
if [ ! -d "scripts" ]; then
    mkdir -p "scripts"
    echo "Created folder: scripts"
fi

# Check and create fastq.list.txt file
if [ ! -f "input/fastq.list.txt" ]; then
    touch "input/fastq.list.txt"
    echo "Created file: input/fastq.list.txt"
fi

echo "Folder structure setup complete."
