#! /bin/sh
# $Id: fwd_check.sh,v 1.42 2008/02/27 17:19:43 lavr Exp $
# Author:   Denis Vakatov (vakatov@ncbi,nlm.nih.gov)
# Modified: Anton Lavrentiev (lavr@ncbi.nlm.nih.gov)
#
# Check for the status of FWDAEMON on the dispatching NCBI servers

delay_sec="$1"
delay_sec=${delay_sec:="10"}

cat <<EOF
http://www.ncbi.nlm.nih.gov/IEB/ToolBox/NETWORK/firewall.html

Checking connections to NCBI Firewall Daemons as of `date -u +'%b %d %Y %R GMT'`:
EOF

{
cat <<EOF
;130.14.25.13	5555	RETIRED
10.10.150.7	5555	INTERNAL
130.14.29.112	5860	RESERVED
130.14.29.112	5861	RESERVED
130.14.29.112	5862	RESERVED
130.14.29.112	5863	RESERVED
130.14.29.112	5864	OK
130.14.29.112	5865	OK
130.14.29.112	5866	RESERVED
130.14.29.112	5867	OK
130.14.29.112	5868	OK
130.14.29.112	5869	OK
130.14.29.112	5870	RESERVED
130.14.29.112	4444	RETIRED
130.14.29.112	4445	RESERVED
130.14.29.112	4446	RESERVED
130.14.29.112	4447	RESERVED
130.14.29.112	4448	OK
130.14.29.112	4449	OK
130.14.29.112	4450	RESERVED
130.14.29.112	4451	OK
130.14.29.112	4452	OK
130.14.29.112	4453	OK
130.14.29.112	4454	RESERVED
EOF
} |
while read x_host x_port x_status ; do
    test "`echo $x_host | grep -c '^[;]'`" != "0"  &&  continue
    if [ "$x_port" -lt "5860" -o "$x_port" -gt "5870" ]; then
        test -z "$HTTP_CAF" -o -n "$HTTP_CAF_EXTERNAL"  &&  continue
    fi
    if [ "$x_status" = "RETIRED"  -o \
         "$x_status" = "RESERVED" ]; then
        echo "${x_host}:${x_port}	$x_status"
        continue
    fi
    test "$x_status" = "READYING"  &&  unset x_status
    ( echo ; sleep $delay_sec ) | telnet $x_host $x_port >/tmp/$$ 2>&1 &
    pid=$!
    trap 'rm -f /tmp/$$; kill $pid >/dev/null 2>&1' 1 2 15
    ( sleep `expr $delay_sec + 2`  &&  kill $pid ) >/dev/null 2>&1 &
    guard=$!
    wait $pid >/dev/null 2>&1
    kill $guard >/dev/null 2>&1
    test -n "$HTTP_CAF_EXTERNAL"  || \
        cp="`tail +4 /tmp/$$ 2>/dev/null | grep -s '^[0-9]\{1,3\}[.][0-9]\{1,3\}[.][0-9]\{1,3\}[.][0-9]\{1,3\}:[0-9]\{1,5\}'`"
    grep -qs 'NCBI Firewall Daemon:  Invalid ticket\.  *Connection closed\.' /tmp/$$ >/dev/null 2>&1
    if   [ $? -eq 0 ]; then
        echo "${x_host}:${x_port}	${x_status:-OKAY}${cp:+	}${cp}"
    elif [ -z "$x_status" ]; then
        echo "${x_host}:${x_port}	READYING"
    else
        echo "${x_host}:${x_port}	FAILED	( telnet $x_host $x_port )"
    fi
    rm -f /tmp/$$
done 2>&1 | grep -v 'Terminated'
