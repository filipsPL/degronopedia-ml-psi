#!/bin/bash

logFile=tests/results.txt
logFileOk=tests/results-ok.txt

echo -n "" > $logFile

for seqfile in $(find seqs/ -type f -name "*.txt" | sort)
do
    echo | tee -a $logFile
    echo "$seqfile" | tee -a $logFile
    cat "$seqfile" | tee -a $logFile

    ./calculate-desc.py --sequence "$seqfile" --type C | tee -a $logFile
    ./calculate-desc.py --sequence "$seqfile" --type NiMetYes | tee -a $logFile
    ./calculate-desc.py --sequence "$seqfile" --type NiMetNo | tee -a $logFile
done



if diff -q $logFile $logFileOk &>/dev/null; then
  echo "Check OK"
  exit 0
else
  echo "Not good!"
  exit 1
fi

