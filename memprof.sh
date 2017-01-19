#!/bin/bash

CMD="genbasis.py 10 1 12"
FNAME="memtest.out"

echo $CMD >> $FNAME
python -m memory_profiler $CMD >> $FNAME

