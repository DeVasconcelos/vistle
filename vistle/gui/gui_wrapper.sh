#! /bin/bash

GDB=""
if [ "$1" = "-g" ]; then
   GDB="gdb --args"
   shift
fi

export DYLD_LIBRARY_PATH="$VISTLE_DYLD_LIBRARY_PATH"

exec $GDB "$(dirname $0)/Vistle.app/Contents/MacOS/Vistle" "$@"
