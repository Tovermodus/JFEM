#!/bin/sh
mvn clean compile exec:java -file='Examples/pom.xml' -Dexec.mainClass='dlm.DLMBenchmark'
