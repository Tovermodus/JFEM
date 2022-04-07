FROM openjdk:17-alpine
RUN apk add maven
RUN apk add git
RUN apk add openblas
RUN apk add openblas-dev
RUN apk add lapack
RUN apk add lapack-dev
RUN mkdir /app
RUN mkdir /app/dlm
RUN echo "hi"
WORKDIR /app
RUN git clone https://github.com/Tovermodus/JFEM.git
WORKDIR /app/JFEM
RUN mvn clean -file=Core/pom.xml
RUN mvn install -file=Core/pom.xml -DskipTests
