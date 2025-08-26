#!/usr/bin/bash

BASEDIR=$(dirname $0)
cd $BASEDIR/frontend
npm start & python ../backend/main.py
