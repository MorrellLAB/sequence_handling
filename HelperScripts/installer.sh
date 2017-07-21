#!/bin/bash

set -e
set -o pipefail

#   This script installs Sickle, Seqqs, and Scythe
#   from GitHub. No GitHub account is required for
#   this installation to occur

#   This also can install other dependencies required
#   by other programs in this repository

#   Create a software directory for these three programs in home directory
#   These paths can be changed to match your hierarchy
#   Note where Seqqs is installed, this is needed in
#   QSub_trim_autoplot.sh to tell where to look for
#   the quality trimming script

usage() {
    echo -e "\
Usage: ./installer.sh <option/> \n\
\n\
Where <option/> is one of 'install', 'bioawk', 'samtools', or 'R' \n\

'./installer.sh install'  installs 'Sickle', 'Scythe' and 'Seqqs' from GitHub for use \n\
with 'Quality_Trimming.sh' script. This script will spit out the path for the \n\
'Seqqs' directory. This needs to be put into 'Quality_Trimming.sh' before use. \n\
\n\
'./installer.sh bioawk' installs 'Bioawk' from GitHub \n\
for use with the 'read_counts.sh' script. \n\
\n\
'./installer.sh samtools' installs 'SAMTools' from GitHub \n\
for use with the 'Read_Mapping.sh', 'SAM_to_BAM.sh', and 'Coverage_Map' scripts. \n\
\n\
'./installer.sh R' installs R from CRAN \n\
for use with the 'Quality_Trimming.sh' and 'Plot_Coverage.sh' scripts. \n\
For those running on MSI, there is built-in loading of the R module. \n\
NOTE: This requires 'wget' to be installed for fetching R from CRAN \n\
\n\
NOTE: Neither GNU Parallel nor the Burrows-Wheeler Aligner (BWA) \n\
are installed using this script. \n\
GNU Parallel is used in every script with capitals in the name. \n\
BWA is used for the 'read_mapping_start.sh' script. \n\
If not using the Minnesota Supercomputing Institute's resources, \n\
you will need to install these manually from their websites. \n\
\n\
BWA:            http://bio-bwa.sourceforge.net/
GNU Parallel:   http://www.gnu.org/software/parallel/
" >&2
    exit 1
}

case "$1" in
    "install" )
        #   Create a software directory in the home directory
        ROOT="${HOME}"
        if [ -d "${ROOT}"/software ]; then
            cd "${ROOT}"/software
            SOFT=`pwd`
        else
            cd "${ROOT}"
            mkdir software
            cd software
            SOFT=`pwd`
        fi
        #   Install Seqqs from the MorrellLab Seqqs GitHub Repository
        echo "Fetching Seqqs from GitHub"
        cd "${SOFT}"
        git clone https://github.com/MorrellLab/seqqs.git
        cd seqqs
        echo "Installing Seqqs"
        make
        echo "Done"
        SEQQS_DIR=`pwd`
        echo export PATH='$PATH':"${SEQQS_DIR}" >> "${HOME}"/.bash_profile
        #   Install Sickle from Vince Buffalo's Sickle GitHub Repository
        echo "Fetching Sickle from GitHub"
        cd "${SOFT}"
        git clone https://github.com/vsbuffalo/sickle.git
        cd sickle
        echo "Installing Sickle"
        make
        echo "Done"
        SICKLE_DIR=`pwd`
        echo export PATH='$PATH':"${SICKLE_DIR}" >> "${HOME}"/.bash_profile
        #   Install Scythe from Vince Buffalo's Scythe GitHub repository
        echo "Fetching Scythe form GitHub"
        cd "${SOFT}"
        git clone https://github.com/vsbuffalo/scythe.git
        cd scythe
        echo "Installing Scythe"
        make all
        echo "Done"
        SCYTHE_DIR=`pwd`
        echo export PATH='$PATH':"${SCYTHE_DIR}" >> "${HOME}"/.bash_profile
        cd "${ROOT}"
        echo "Seqqs directory is ${SEQQS_DIR}"
        sleep 2
        echo 'This needs to be written in "QSub_trim_autoplot.sh"'
        sleep 3
        echo "To finish installing, please run the following command:"
        echo "source ~/.bash_profile"
        ;;
    "bioawk" )
        #   Create a software directory in the home directory
        ROOT=`pwd`
        if [ -d "${ROOT}"/software ]; then
            cd "${ROOT}"/software
            SOFT=`pwd`
        else
            cd "${ROOT}"
            mkdir software
            cd software
            SOFT=`pwd`
        fi
        #   Check to see if Bioawk is already installed
        if `command -v bioawk > /dev/null 2> /dev/null`
        then
            echo "Bioawk is installed"
        else
            #   If not, download Bioawk from GitHub and install
            cd "${SOFT}"
            git clone https://github.com/lh3/bioawk.git
            cd bioawk
            make
            BIOAWK_DIR=`pwd`
            echo export PATH='$PATH':"${BIOAWK_DIR}" >> "${HOME}"/.bash_profile
        fi
        cd "${ROOT}"
        echo "To finish installing, please run the following command:"
        echo "source ~/.bash_profile"
        ;;
    "samtools" )
        #   Create a software directory in the home directory
        ROOT=`pwd`
        if [ -d "${ROOT}"/software ]; then
            cd "${ROOT}"/software
            SOFT=`pwd`
        else
            cd "${ROOT}"
            mkdir software
            cd software
            SOFT=`pwd`
        fi
        #   Check to see if SAMTools is already installed
        if `command -v samtools > /dev/null 2> /dev/null`
        then
            echo "sammtools is installed"
        else
            #   If not, download SAMTools from GitHub and install
            cd "${SOFT}"
            git clone https://github.com/samtools/htslib.git
            cd htslib
            autoconf
            ./configure --prefix=`pwd`
            make
            make install
            HTSLIB_DIR=`pwd`
            echo export PATH='$PATH':"${HTSLIB_DIR}" >> "${HOME}"/.bash_profile
            cd "${SOFT}"
            git clone https://github.com/samtools/samtools.git
            cd samtools
            make
            SAMTOOLS_DIR=`pwd`
            echo export PATH='$PATH':"${SAMTOOLS_DIR}" >> "${HOME}"/.bash_profile
        fi
        cd "${ROOT}"
        echo "To finish installing, please run the following command:"
        echo "source ~/.bash_profile"
        ;;
    "R" )
        #   Create a software directory in the home directory
        ROOT=`pwd`
        if [ -d "${ROOT}"/software ]; then
            cd "${ROOT}"/software
            SOFT=`pwd`
        else
            cd "${ROOT}"
            mkdir software
            cd software
            SOFT=`pwd`
        fi
        #   Check to see if R is already installed
        if `command -v Rscript > /dev/null 2> /dev/null`
        then
            echo "R is installed"
        else
            #   If no R, check to see if wget is installed
            if `command -v wget > /dev/null 2> /dev/null`
            then
                #   If no R, but yes wget, download R from CRAN and install
                cd "${SOFT}"
                wget --no-directories --progress=bar -r -A.tar.gz http://cran.r-project.org/
                rm robots.txt
                tar -xvzf *.tar.gz
                cd R*
                ./configure --prefix=`pwd`
                make
                cd bin
                R_DIR=`pwd`
                echo export PATH='$PATH':"${R_DIR}" >> "${HOME}"/.bash_profile
            else
                #   If neither R nor wget are installed, exit out
                echo "Please install wget to your system for to install R using this script"
                exit 1
            fi
        fi
        cd "${ROOT}"
        echo "To finish installing, please run the following command:"
        echo "source ~/.bash_profile"
        ;;
    * )
        usage
        ;;
esac
