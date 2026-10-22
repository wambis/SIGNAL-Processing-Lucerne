#!/bin/bash

PID=628042

echo "Waiting for PID $PID to finish..."
while kill -0 $PID 2>/dev/null; do sleep 5; done

echo "Starting next job..."
nohup Plot_coherence_SEVERAL-ALL-2D.py > log_2D.out 2>&1 &
