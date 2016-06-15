#/bin/sh

SRC_FILES = src/stiffness.f90 \
	src/stress.f90

default:
	f2py --fcompiler=gnu95 -c ${SRC_FILES} -m lib
