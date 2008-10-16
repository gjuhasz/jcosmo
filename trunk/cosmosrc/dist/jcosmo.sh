#!/bin/sh
cd `dirname $0`
if [ -n "$JAVA_HOME" ]; then
  $JAVA_HOME/bin/java -jar ./cosmo.jar $*
else
  java -jar ./jcosmo.jar $*
fi
cd $OLDPWD
