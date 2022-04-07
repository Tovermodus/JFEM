FROM openjdk:17-alpine
RUN apk add maven
RUN apk add git
RUN apk add openblas
RUN apk add openblas-dev
RUN apk add lapack
RUN apk add lapack-dev
RUN rm /etc/mavenrc
RUN echo 'M2_HOME="$m2_home"' > /etc/mavenrc
RUN mkdir /app
RUN mkdir /app/dlm
RUN echo "hi"
ENV MAVEN_OPTS="-Xmx8000m"
WORKDIR /app
RUN git clone https://github.com/Tovermodus/JFEM.git
WORKDIR /app/JFEM/JSparse
RUN /bin/sh build.sh
RUN /bin/sh run.sh
RUN cp out/artifacts/jSparse.jar ../jSparse.jar
WORKDIR /app/JFEM
RUN mvn clean -file=Core/pom.xml
RUN mvn install -file=Core/pom.xml -DskipTests
