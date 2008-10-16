#!/bin/sh
cd `dirname $0`
if [ -n "$JAVA_HOME" ]; then
  $JAVA_HOME/bin/java -jar ./cosmo.jar $*
else
  java -jar ./cosmo.jar $*
fi
cd $OLDPWD
