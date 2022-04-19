# Dockerfile sreenaths 18/04/2022
FROM alpine
RUN apk update
# get the make and linux-headers to compile gdb
RUN apk add --no-cache make
RUN apk add --no-cache linux-headers
RUN apk add --no-cache texinfo
RUN apk add --no-cache gcc
RUN apk add --no-cache g++
RUN apk add --no-cache gfortran
# install gdb
RUN apk add --no-cache gdb
# install python and fortran language-server
RUN apk add --no-cache python3 py3-pip
RUN pip install -U setuptools
RUN pip install fortran-language-server
RUN apk add --no-cache openssh
RUN ssh-keygen -A
# use foo user for ssh login
RUN adduser -D naoe5
RUN echo naoe5:osaka | chpasswd
CMD ["/usr/sbin/sshd","-D"]