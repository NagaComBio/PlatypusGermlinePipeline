#!/bin/bash
#
# Program to call the Platypus Pipeline
#  by Nagarajan Paramasivam
#  Nov 2016
# #############################################################################
# 
## For cancer samples: control and tumor BAMs for PID_1
#
#  $ sh run_platypusPipeline.sh -c CONFIG_FILE PID_1
# -----------------------------------------------------------------------------
## For TRIO samples or Familial: Normal/control sample for all the PIDs in a 
## family
#### First pid after the config file path will be selected as representative
#### pid. The ouput VCF will be stored in rep.PID. If separate option is 
#### provided in the config file, then softlink will be created in other pids
#### The separate option should be switched on, if you are planning for the 
#### downstream Germline Analysis
#  $ sh run_platypusPipeline.sh -c CONFIG_FILE PID_1 PID_2 PID_3
##############################################################################

#### Usage function
function Exit_Usage {
  printf "Usage: %s -c CONFIG_FILE PID_1 PID_2 ...\n" $(basename $0) >&2
  exit 2
}

### Argument loop
# If no argument is provide exit with usage message
if [ $# -eq 0 ]
then
  Exit_Usage
fi

while getopts ':c:' OPTION
do
  case $OPTION in
    c) cflag=1
       CONFIG_FILE="$OPTARG"
    ;;
    *) Exit_Usage
    ;;
  esac
done


# Everything after the flagged arguments are PIDs
shift $(($OPTIND - 1))
ALL_PIDS=$*

# If there are no pids then exit with usage 
if [[ -z $ALL_PIDS ]]
then
  Exit_Usage
fi

# Selecting first PID as representative pid
joined_PIDS=$(echo $ALL_PIDS | tr " " "#")

tempPIDS=($(echo $joined_PIDS | tr "#" " "))
PID=${tempPIDS[0]}  # reprenstative PID

# Reading in Configure and copying it to representative pid subfolders
source $CONFIG_FILE

if [[ ! -d $ANALYSIS_FOLDER/$PID/$SUBDIR ]]
then 
  mkdir -p -m 770 $ANALYSIS_FOLDER/$PID/$SUBDIR
fi

cp $CONFIG_FILE $ANALYSIS_FOLDER/$PID/$SUBDIR/

###############################################################################
# Running Platypus and Splitting the raw file into SNVs and Indels
#
platypusRun_depend=''
echo -e "$PID:"

if [[ $RUN_PLATYPUS == 1 ]]
  then    
    platypusRun=`bsub -q verylong -J platypusRun_${PID} -n 10 -R "span[ptile=10]" -R "rusage[mem=10480]" -env "all,CONFIG_FILE=$CONFIG_FILE,PID=$PID,PIDs=$joined_PIDS" sh $PIPELINE_DIR/platypusIndel.sh | cut -d'<' -f2 | cut -d'>' -f1`
    platypusRun_depend="-w \"done($platypusRun)\""
    echo -e "    Platypus submitted    : $platypusRun"
  else
  	platypusRun=""
fi

sleep 2
###############################################################################
# IndelAnnotation
#
if [[ "$RUN_INDELANNOTATION" == 1 ]]
  then    
    indelAnnotator=`bsub -q verylong -J indelAnnotation_${PID} -R "rusage[mem=5480]" $platypusRun_depend -env "all,CONFIG_FILE=$CONFIG_FILE,PID=${PID}" sh $PIPELINE_DIR/indelAnnotation.sh | cut -d'<' -f2 | cut -d'>' -f1`    
    echo -e "    Indel annotation      : $indelAnnotator"
fi


################################################################################
# SNV Annotation
#
if [[ "$RUN_SNVANNOTATION" == 1 ]]
  then            
    snvAnnotator=`bsub -q verylong -J snvAnnotation_${PID} -R "rusage[mem=5480]" $platypusRun_depend -env "all,CONFIG_FILE=$CONFIG_FILE,PID=${PID}" sh $PIPELINE_DIR/snvAnnotation.sh | cut -d'<' -f2 | cut -d'>' -f1`
    echo -e "    SNV annotation        : $snvAnnotator"
fi


###############################################################################
# Creating softlinks to all the PIDs in multi-calling
#
if [[ "$RUN_SEPARATE_PIDs" == 1 ]]
  then    
    separatorsPIDs=`bsub -q short -J separatePIDs_${PID} -w "done($indelAnnotator) && done($snvAnnotator)" -env "all,CONFIG_FILE=$CONFIG_FILE,PIDs=$joined_PIDS,PID=$PID" sh $PIPELINE_DIR/separate_PIDs.sg | cut -d'<' -f2 | cut -d'>' -f1`    
    echo -e "    Separate PIDs   : $separatorsPIDs"
fi