#!/bin/bash

EVENTS=1000000

for Ex in 0.0 0.1 0.5 2.5
do
  ./exec_UNISim-tool.exe -events ${EVENTS} -o TheOutput_Ex${Ex}.root -physics SequentialDecay -reaction ./reactions/TheReactionFile_Ex${Ex}.dat
done
