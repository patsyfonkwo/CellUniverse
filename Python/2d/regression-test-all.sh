#!/bin/bash
# Functions
die(){ (echo "$USAGE"; echo "FATAL ERROR: $@")>&2; exit 1; }
warn(){ (echo "WARNING: $@")>&2; }
not(){ if eval "$@"; then return 1; else return 0; fi; }
newlines(){ awk '{for(i=1; i<=NF;i++)print $i}' "$@"; }
parse(){ awk "BEGIN{print $@}" </dev/null; }
cpus() {
    TMP=/tmp/cpus.$$
    trap "/bin/rm -f $TMP; exit" 0 1 2 3 15

    # Most Linux machines:
    lscpu >$TMP 2>/dev/null && awk '/^CPU[(s)]*:/{cpus=$NF}END{if(cpus)print cpus; else exit 1}' $TMP && return

    # MacOS:
    ([ `arch` = Darwin -o `uname` = Darwin ] || uname -a | grep Darwin >/dev/null) && sysctl -n hw.ncpu && return

    # Cygwin:
    case `arch` in
    CYGWIN*) grep -c '^processor[ 	]*:' /proc/cpuinfo; return ;;
    *) if [ -d /dev/cpu -a ! -f /dev/cpu/microcode ]; then
	ls -F /dev/cpu | fgrep -c
	return
       fi
	;;
    esac

    # Oops
    echo "couldn't figure out number of CPUs" >&2; exit 1
}

# generally useful Variables
NL='
'
TAB='	'

# Temporary Filename + Directory (both, you can use either, note they'll have different random stuff in the XXXXXX part)
TMP=`mktemp /tmp/$BASENAME.XXXXXX`
TMPDIR=`mktemp -d /tmp/$BASENAME.XXXXXX`
trap "/bin/rm -rf $TMP $TMPDIR; exit" 0 1 2 3 15 # call trap "" N to remove the trap for signal N

while [ $# -gt 0 ]; do
    case "$1" in
    -use-git-at)
	if [ -f git-at ] && [ `wc -l < git-at` -eq 2 -a `git log -1 --format=%at` -eq `tail -1 git-at` ]; then
	    echo -n "Repo unchanged; returning same status code as "
	    tail -1 git-at | xargs -I{} date -d @{} +%Y-%m-%d-%H:%M:%S
	    exit `head -1 git-at`
	fi
	;;
    -make) echo "Nothing to make" >&2;;
    *) break;;
    esac
    shift
done

USAGE="USAGE: $0 [ list of tests to run, defaults to regression-tests/*/*.sh ]"

PATH=`pwd`:`pwd`/scripts:$PATH
export PATH

CORES=${CORES:=`cpus 2>/dev/null | awk '{c2=int($1/2); if(c2>0)print c2; else print 1}'`}
[ "$CORES" -gt 0 ] || die "can't figure out how many cores this machine has"
echo "Using $CORES cores for regression tests"
export EXE CORES

NUM_FAILS=0
STDBUF=''
if which stdbuf >/dev/null; then
    STDBUF='stdbuf -oL -eL'
fi
if [ $# -eq 0 ]; then
    set regression-tests/*/*.sh
fi
for r
do
    REG_DIR=`dirname "$r"`
    NEW_FAILS=0
    export REG_DIR
    echo --- running test $r ---
    if eval time $STDBUF "$r"; then # force output and error to be line buffered
	:
    else
	NEW_FAILS=$?
	(( NUM_FAILS+=$NEW_FAILS ))
    fi
    echo --- test $r incurred $NEW_FAILS failures, cumulative failures is $NUM_FAILS ---
done
echo Total number of failures: $NUM_FAILS
(echo $NUM_FAILS; git log -1 --format=%at) > git-at
exit $NUM_FAILS
