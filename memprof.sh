#!/bin/bash

CMD="genbasis.py 10 1 15"
FNAME="memtest.out"

echo $CMD >> $FNAME
python -m memory_profiler $CMD >> $FNAME

