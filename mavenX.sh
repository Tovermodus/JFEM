#!/bin/sh
export MAVEN_OPTS="-Xmx30000m"
mvn clean compile exec:java -file='Examples/pom.xml' -Dexec.mainClass='dlm.DLMBenchmark'
