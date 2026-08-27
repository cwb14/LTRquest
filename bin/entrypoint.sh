#!/bin/sh
# Entry point for the LTRquest container image.
#
# Apptainer turns a Docker ENTRYPOINT into the SIF's runscript, which is what
# decides whether `./ltrquest.sif --genome x.fa` means anything. An image with
# only a CMD gets a runscript that *replaces* the command with whatever the user
# typed, so 1.0.1 answered
#
#     ./ltrquest.sif --help
#     FATAL: "--help": executable file not found in $PATH
#
# Setting an ENTRYPOINT makes Apptainer prepend it to the user's arguments
# instead, and this script decides what those arguments mean:
#
#     ./ltrquest.sif --genome x.fa       -> ltrquest --genome x.fa
#     ./ltrquest.sif ltrquest-gff3 -h    -> ltrquest-gff3 -h
#     docker run IMG ltrquest --help     -> ltrquest --help      (as documented)
#     docker run IMG /bin/bash -c '...'  -> /bin/bash -c '...'   (Nextflow)
#     apptainer exec IMG ltrquest ...    -> never reaches here; exec skips it
#
# The rule: a first argument that names a runnable command is a command;
# anything else is an ltrquest argument. That is strictly more permissive than
# the empty ENTRYPOINT it replaces, so no invocation that worked before stops
# working.

set -eu

case "${1-}" in
    ""|-*)
        # No arguments, or the first one is a flag: both belong to ltrquest.
        # Matched before `command -v` because `command -v --help` is not
        # portable across shells.
        ;;
    *)
        if command -v "$1" > /dev/null 2>&1; then
            exec "$@"
        fi
        ;;
esac

exec ltrquest "$@"
