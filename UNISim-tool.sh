#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
export UNISimSrc=${DIR}/

export LD_LIBRARY_PATH=$UNISimSrc/lib:$LD_LIBRARY_PATH
