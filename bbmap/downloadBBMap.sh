#!/bin/bash

wd=`pwd`
ad=${wd}"/bbmap"
cd ${ad}

# Script to download and install BBMap suite
# BBMap is a collection of bioinformatics tools including clumpify
# for deduplication

set -e

echo "================================================"
echo "BBMap Suite Installation Script"
echo "================================================"

# Version to download
BBMAP_VERSION="39.46"
BBMAP_URL="https://sourceforge.net/projects/bbmap/files/BBMap_${BBMAP_VERSION}.tar.gz/download"
BBMAP_TAR="BBMap_${BBMAP_VERSION}.tar.gz"

echo ""
echo "Download directory: ${ad}"
echo "BBMap version: ${BBMAP_VERSION}"
echo ""

# Check if Java is installed
if ! command -v java &> /dev/null; then
    echo "WARNING: Java not found. BBMap requires Java 8 or higher."
    echo "Please install Java before running BBMap tools."
    echo ""
fi

# Download BBMap
echo "Downloading BBMap ${BBMAP_VERSION}..."
wget -O "${BBMAP_TAR}" "${BBMAP_URL}"

if [ ! -f "${BBMAP_TAR}" ]; then
    echo "ERROR: Download failed. File not found: ${BBMAP_TAR}"
    exit 1
fi

# Extract the archive
echo "Extracting BBMap..."
tar -xzf "${BBMAP_TAR}"

if [ ! -d "bbmap" ]; then
    echo "ERROR: Extraction failed. Directory 'bbmap' not found."
    exit 1
fi

# Make shell scripts executable
echo "Setting permissions..."
chmod +x bbmap/*.sh 2>/dev/null || true

# Clean up
echo "Cleaning up..."
mv bbmap/* .
rm -rf bbmap
rm -f "${BBMAP_TAR}"

# Verify installation
echo ""
echo "================================================"
echo "Installation Complete!"
echo "================================================"
echo ""
echo "BBMap installed at: ${ad}"
echo ""

# Check key tools
if [ -f "clumpify.sh" ]; then
    echo "✓ clumpify.sh found"
else
    echo "✗ clumpify.sh not found"
fi

echo ""
echo "Documentation is available in: ${ad}/docs/"
echo ""
echo "To test the installation, run:"
echo "  ${ad}/clumpify.sh --version"
echo ""
echo "Note: Make sure Java 8 or higher is installed on your system."
echo "================================================"

cd ${wd}

exit
