FROM archlinux
RUN pacman -Syu --noconfirm
RUN pacman -S jdk11-openjdk gcc maven git openblas --noconfirm
RUN pacman -S lapack --noconfirm
#RUN rm /etc/mavenrc
#RUN echo 'M2_HOME="$m2_home"' > /etc/mavenrc
RUN mkdir /app
RUN mkdir /app/dlm
ENV MAVEN_OPTS="-Xmx8000m"
WORKDIR /app
COPY . /app/JFEM
#RUN git clone https://github.com/Tovermodus/JFEM.git
WORKDIR /app/JFEM/JSparse
#RUN echo "his"
RUN ./build.sh
RUN ./run.sh
RUN cp /app/JFEM/JSparse/out/artifacts/jSparse/jSparse.jar ../jSparse.jar
WORKDIR /app/JFEM
RUN mvn clean -file=Core/pom.xml
RUN mvn install -file=Core/pom.xml -DskipTests
