#!/bin/bash

# Copyright 2013-2015, Derrick Wood, Jennifer Lu <jlu26@jhmi.edu>
#
# This file is part of the Kraken taxonomic sequence classification system.
#
# Kraken is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Kraken is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Kraken.  If not, see <http://www.gnu.org/licenses/>.

# Check that jellyfish is executable and is proper version
# Also check for KMC tools if KRAKEN_USE_KMC is set
# Designed to be called by kraken-build

set -u  # Protect against uninitialized vars.
set -e  # Stop on error
set -o pipefail  # Stop on failures in non-final pipeline commands

# Check if we should use KMC instead of Jellyfish
if [ -n "${KRAKEN_USE_KMC:-}" ]
then
  echo "KMC mode enabled, checking for KMC tools..."
  
  # Check for kmc
  if ! command -v kmc &> /dev/null; then
    echo "ERROR: kmc command not found in PATH"
    echo "Please install KMC3 and ensure it's in your PATH"
    exit 1
  fi
  
  # Check for kmc_tools
  if ! command -v kmc_tools &> /dev/null; then
    echo "ERROR: kmc_tools command not found in PATH"
    echo "Please install KMC3 and ensure kmc_tools is in your PATH"
    exit 1
  fi
  
  # Check for kmc_to_jellyfish
  if [ -f "$(dirname "$0")/../src/kmc_to_jellyfish" ]; then
    echo "Found kmc_to_jellyfish in src directory"
  else
    echo "ERROR: kmc_to_jellyfish tool not found in src directory"
    echo "Please build the kmc_to_jellyfish tool with 'make -C src'"
    exit 1
  fi
  
  echo "Found KMC tools: kmc, kmc_tools, kmc_to_jellyfish"
else
  # Check for Jellyfish
  JELLYFISH_VERSION=$(jellyfish --version | awk '{print $2}')
  if [[ $JELLYFISH_VERSION =~ ^1\. ]]
  then
    echo "Found jellyfish v$JELLYFISH_VERSION"
  else
    echo "Found jellyfish v$JELLYFISH_VERSION"
    echo "Kraken requires jellyfish version 1"
    exit 1
  fi
fi
