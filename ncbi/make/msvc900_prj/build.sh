#! /bin/sh
# $Id: build.sh,v 1.1 2009/01/23 19:48:09 ivanov Exp $
# Author:  Vladimir Ivanov (ivanov@ncbi.nlm.nih.gov)
#
# Build C Toolkit on MSVC9 32/64.


########### Arguments

script="$0"
cfgs="${1:-Debug Release DebugDLL ReleaseDLL}"
arch="$2"
target="${3:-all_ncbi}"


########### Global variables

timer="date +'%H:%M'"
out=".build.$$"


########## Functions

error()
{
  echo "[`basename $script`] ERROR:  $1"
  exit 1

}

generate_msvc9_error_check_file() {
  cat <<-EOF >$1
	/(^| : |^The source )([fatal error]* [CDULNKPRJVT]*[0-9]*: |The .* are both configured to produce |Error executing )/ {
	  print \$0
	  exit
	}
	EOF
}


########## Main

# Get build dir
build_dir=`dirname $script`
build_dir=`(cd "$build_dir"; pwd)`

if [ ! -d $build_dir ] ; then
  error "Build directory $build_dir not found"
  exit 1
fi
cd $build_dir


# Generate errors check script

check_awk=$build_dir/build_check.awk
generate_msvc9_error_check_file $check_awk


# Build

for cfg in $cfgs ; do
   start=`eval $timer`
   echo Start time: $start
   echo "INFO: Building \"$cfg\""
   $build_dir/build_exec.bat "ncbi.sln" build "$arch" "$cfg" "$target" $out
   status=$?
   cat $out
   echo "Build time: $start - `eval $timer`"
   if [ $status -ne 0 ] ; then
     # Check on errors
     failed="1"
     grep '^ *Build: .* succeeded, .* failed' /tmp/build.$$ >/dev/null 2>&1  && \
       awk -f $check_awk $out >$out.res 2>/dev/null  &&  test ! -s $out.res  &&  failed="0"
     rm -f $out.res >/dev/null 2>&1
     rm -f $out     >/dev/null 2>&1
     if [ "$failed" = "1" ]; then
       exit 4
     fi
   fi
   rm -f $out >/dev/null 2>&1
done


exit 0
