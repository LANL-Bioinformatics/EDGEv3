#!/bin/bash

export RATT_HOME=$CONDA_PREFIX/opt/RATT #remove when Dockerizeing

$RATT_HOME/start.ratt.sh $1 $2 $3 $4