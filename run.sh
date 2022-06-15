#!/usr/bin/env bash

# out of source builds
cmake -S ${GITPOD_REPO_ROOT}/flow-test-dumux -B ${GITPOD_REPO_ROOT}/build
cmake --build ${GITPOD_REPO_ROOT}/build
# CASES=(tracer_transport hello_world freeflow cavity shallow_water rotation infiltration porenetwork) # ${CASES[@]}
tar -c -I 'zstd -19 -T0' -f ${GITPOD_REPO_ROOT}/build_$(date -u +"%Y-%m-%d-%H-%M-%S" --date='5 hours ago').tar.zst ${GITPOD_REPO_ROOT}/build/{tracer_transport,hello_world,freeflow,cavity,shallow_water,rotation,infiltration,porenetwork}