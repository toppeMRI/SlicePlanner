
## Install jdk 12
 - Downloaded JDK 12 from http://jdk.java.net/12/
 - openjdk-12.0.2_linux-x64_bin.tar.gz
 - untar

```
sudo mv jdk-12.0.2/ /usr/lib/jvm/
cd /etc/alternatives/
sudo rm -f java
sudo ln -s /usr/lib/jvm/jdk-12.0.2/bin/java java
sudo rm -f javac
sudo ln -s /usr/lib/jvm/jdk-12.0.2/bin/javac javac
```

In .bashrc
```
export JAVA_HOME=/usr/lib/jvm/jdk-12.0.2
```

## Install JavaFX 12

Getting started with JavaFX: https://openjfx.io/openjfx-docs/

 - Downloaded JavaFX runtime from https://gluonhq.com/products/javafx/
 - openjfx-12.0.2_linux-x64_bin-sdk.zip
 - unzip
```
sudo mv javafx-sdk-12.0.2/ /usr/lib/
```

In .bashrc
```
export PATH_TO_FX=/usr/lib/javafx-sdk-12.0.2/lib
```

## Install HDF5

 - Download sis-jhdf5-19.04.0.zip
 - unzip


## Compile and run 

Run javac and java with the appropriate classpath, modules, etc. See ../sandbox/doit.

```
./doit
```

