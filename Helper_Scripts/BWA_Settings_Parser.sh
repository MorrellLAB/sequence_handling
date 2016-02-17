#!/bin/bash

POSITIONALS='' # Create a string of positional arguments for BWA meme
if [[ "${PAIRED}" == true ]]; then POSITIONALS="${POSITIONALS}"'-P '; fi # Add paired to the positionals
if [[ "${INTERLEAVED}" == true ]]; then POSITIONALS="${POSITIONALS}"'-p '; fi # Add interleaved to the positionals
if [[ "${SECONDARY}" == true ]]; then POSITIONALS='-a'; else SECONDARY=''; fi # Add secondary to the positionals
if [[ "${APPEND}" == true ]]; then POSITIONALS="${POSITIONALS}"'-C '; fi # Add append to the positionals
if [[ "${HARD}" == true ]]; then POSITIONALS="${POSITIONALS}"'-H '; fi # Add hard to the positionals
if [[ "${SPLIT}" == true ]]; then POSITIONALS="${POSITIONALS}"'-M '; fi # Add split to the positionals

if [[ "${VERBOSITY}" == 'disabled' ]]; then VERBOSITY=0; elif [[ "${VERBOSITY}" == 'errors' ]]; then VERBOSITY=1; elif [[ "${VERBOSITY}" == 'warnings' ]]; then VERBOSITY=2; elif [[ "${VERBOSITY}" == 'all' ]]; then VERBOSITY=3; elif [[ "${VERBOSITY}" == 'debug' ]]; then VERBOSITY=4; else echo "Failed to recognize verbosity level, exiting..."; exit 1; fi # Set the verbosity level

MEM_SETTINGS=`echo "-t ${THREADS} -k ${SEED} -w ${WIDTH} -d ${DROPOFF} -r ${RE_SEED} -A ${MATCH} -B ${MISMATCH} -O ${GAP} -E ${EXTENSION} -L ${CLIP} -U ${UNPAIRED} -T ${THRESHOLD} -v ${VERBOSITY} ${POSITIONALS}"` # Store the settings
