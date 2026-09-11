#!/bin/bash

set -ex

export ASTERXSPACE="$PWD"
export WORKSPACE="$PWD/../workspace"
mkdir -p "$WORKSPACE"
cd "$WORKSPACE"

cd Cactus

# SimFactory commit 2718e43 added db.ini with the same [db] section already
# present in db1.hpc.lsu.edu.ini. Its loader rejects duplicate machine names
# before honoring --machine, so this otherwise unrelated entry breaks CI.
# Disable only that exact conflict and leave the workaround dormant once
# SimFactory fixes either machine definition.
if grep -Fxq '[db]' ./simfactory/mdb/machines/db.ini 2>/dev/null && \
        grep -Fxq '[db]' ./simfactory/mdb/machines/db1.hpc.lsu.edu.ini 2>/dev/null; then
    mv ./simfactory/mdb/machines/db.ini \
        ./simfactory/mdb/machines/db.ini.disabled
fi

# Set up SimFactory
cp "$ASTERXSPACE/scripts/actions-$ACCELERATOR-$REAL_PRECISION.cfg" ./simfactory/mdb/optionlists
cp "$ASTERXSPACE/scripts/actions-$ACCELERATOR-$REAL_PRECISION.ini" ./simfactory/mdb/machines
cp "$ASTERXSPACE/scripts/actions-$ACCELERATOR-$REAL_PRECISION.run" ./simfactory/mdb/runscripts
cp "$ASTERXSPACE/scripts/actions-$ACCELERATOR-$REAL_PRECISION.sub" ./simfactory/mdb/submitscripts
cp "$ASTERXSPACE/scripts/asterx.th" .
cp "$ASTERXSPACE/scripts/defs.local.ini" ./simfactory/etc

# Compiler cache: cuts warm CI rebuilds dramatically. The .ccache dir is
# persisted across runs by actions/cache (see .github/workflows/ci.yml).
# No-op if ccache is not installed. We prepend "ccache " to CC and CXX in the
# copied option list (matches ../AnalysisX/scripts/build.sh).
if command -v ccache >/dev/null 2>&1; then
    export CCACHE_DIR="${CCACHE_DIR:-$ASTERXSPACE/.ccache}"
    ccache --max-size=2G
    ccache --zero-stats || true
    sed -i -e 's/^CC = /CC = ccache /' -e 's/^CXX = /CXX = ccache /' \
        "./simfactory/mdb/optionlists/actions-$ACCELERATOR-$REAL_PRECISION.cfg"
fi

# For Formaline
git config --global user.email "asterx@einsteintoolkit.org"
git config --global user.name "Github Actions"

# Build
# The build log needs to be stored for later.
time ./simfactory/bin/sim --machine="actions-$ACCELERATOR-$REAL_PRECISION" build -j $(nproc) sim 2>&1 |
    tee build.log

# Check whether the executable exists and is executable
test -x exe/cactus_sim

# Report cache effectiveness (warm runs should show a high hit rate).
command -v ccache >/dev/null 2>&1 && ccache --show-stats || true
