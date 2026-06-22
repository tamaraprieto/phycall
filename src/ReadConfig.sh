#!/bin/sh
# Read parameters from a configuration file
#
# Usage:  source ReadConfig.sh /path/to/Config.<PATIENT>.txt
# Help:   ReadConfig.sh --help
#
# This script must be sourced (not executed) so the exported variables are
# available to the calling script. It must live in a directory on $PATH, or be
# sourced with an explicit path (e.g. source "$SCRIPTDIR/ReadConfig.sh" "$1").

if [ "$#" -ne 1 ]
then
        echo "You must specify absolute path to a configuration file as argument. More info about how to create this file: ReadConfig.sh --help"
        echo ""
        echo "--------------------------------------------------------"
        exit
elif [ "$1" = "--help" ]
then
        echo ""
        echo "------------------------------------------------------------------------------------------------------------------------"
        echo ""
        echo "INSTRUCTIONS IN HOW TO CREATE A CONFIG FILE"
        echo ""
        echo "The config file must contain following arguments and directories have to be absolute (environmental variables are not allowed):"
        echo "id | id_list: File containing read group ID samples."
        echo "s | sample_list: File containing one sample name per row. This file must be located in the oridir directory."
        echo "b | bulk_list: File containing the name of the bulk samples."
        echo "c | control_list: File containing the name of the control samples."
        echo "ori | original_directory: Absolute path of the folder containing raw data and sample lists."
        echo "work | working_directory: Absolute path to the working directory."
        echo "res | resources_directory: Absolute path of the reference and other resources directory."
        echo "script | scripts_directory: Absolute path to the scripts directory."
        echo "ref | reference_name: Name of the reference genome."
        echo "targ | target_capture: Name of the file containing sequences selected in a WES experiment or any other kind of targeted sequencing experiment."
        echo "ad | library_adapters: Tab delimited file containing the name of library in the first column and the sequence of the adapters in the other columns."
        echo "lib | library: Name of the sequencing library."
        echo "pl | platform: Name of the sequencing technology among ILLUMINA,SLX,SOLEXA,SOLID,454,LS454,COMPLETE,PACBIO,IONTORRENT,CAPILLARY or HELICOS. If other specify UNKNOWN."
        echo "wga_lib | whole_genome_amplification_library: Name of the WGA kit among AMPLI-1,PICOPLEX or MALBAC."
        echo "e | email: email address to notify failing."
        echo "p | queue: Name of the queue the jobs will be sent to."
        echo "g | gender: Gender of the sample."
        echo ""
        echo "Example of lines in config file:"
        echo ""
        echo "sample_list=mysamples.txt"
        echo "c=mycontrols.txt"
        echo "resources_directory=/mnt/lustre/scratch/home/uvi/be/phylocancer/RESOURCES/"
        echo ""
        echo "These variables will correspond to:"
        echo '$IDLIST,$SAMPLELIST,$BULKLIST,$CONTROL,$ORIDIR,$WORKDIR,$RESDIR,$SCRIPTDIR,$REF,$TARGET,$LIBRARY,$WGA_LIBRARY,$PLATFORM,$ADAPTERS,$QUEUE,$EMAIL,$GENDER'
        echo "-----------------------------------------------------------------------------------------------------------------------"
        echo ""
        exit
elif [ ! -e "$1" ]
then
        echo "File $1 does not exist"
        exit
else
        while read -r i; do
                # Skip blank lines and comments
                case $i in
                ''|\#*) continue ;;
                esac
                case $i in
                id=*|id_list=*)                              export IDLIST="${i#*=}" ;;
                s=*|sample_list=*)                           export SAMPLELIST="${i#*=}" ;;
                b=*|bulk_list=*)                             export BULKLIST="${i#*=}" ;;
                c=*|control_list=*)                          export CONTROL="${i#*=}" ;;
                ori=*|original_directory=*)                  export ORIDIR="${i#*=}" ;;
                work=*|working_directory=*)                  export WORKDIR="${i#*=}" ;;
                res=*|resources_directory=*)                 export RESDIR="${i#*=}" ;;
                script=*|scripts_directory=*)                export SCRIPTDIR="${i#*=}" ;;
                ref=*|reference_name=*)                       export REF="${i#*=}" ;;
                targ=*|target_capture=*)                     export TARGET="${i#*=}" ;;
                ad=*|library_adapters=*)                     export ADAPTERS="${i#*=}" ;;
                lib=*|library=*)                             export LIBRARY="${i#*=}" ;;
                pl=*|platform=*)                             export PLATFORM="${i#*=}" ;;
                wga_lib=*|whole_genome_amplification_library=*) export WGA_LIBRARY="${i#*=}" ;;
                e=*|email=*)                                 export EMAIL="${i#*=}" ;;
                p=*|queue=*)                                 export QUEUE="${i#*=}" ;;
                g=*|gender=*)                                export GENDER="${i#*=}" ;;
                *)
                        echo "$i is not a valid parameter"
                        exit
                        ;;
                esac
        done < "$1"
fi
