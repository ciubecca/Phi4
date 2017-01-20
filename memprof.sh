#!/bin/bash

CMD="eigs.py 10 1 13"
FNAME="profmem.out"

echo $CMD >> $FNAME
python -m memory_profiler $CMD >> $FNAME

