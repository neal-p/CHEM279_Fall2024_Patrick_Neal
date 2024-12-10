FROM ubuntu:22.04 AS build

RUN apt-get update -qq && apt-get install -qq \
    git \
    cmake \
    gfortran \
    libopenblas-dev \
    python3 \
    python3-pip 

RUN pip install numpy

RUN git clone --depth=1 https://github.com/openmopac/mopac.git

WORKDIR /mopac/build
RUN cmake .. && make


#############################################################################


FROM ubuntu:22.04 AS install

RUN apt-get update -qq && apt-get install -qq \
    gfortran \
    libopenblas-dev \
    python3.11 \
    python3-pip \
    openbabel 

RUN pip install --no-input numpy pandas scikit-learn torch tqdm matplotlib scipy

WORKDIR /mopac/build
COPY --from=build /mopac/build/*mod /mopac/build/
COPY --from=build /mopac/build/*.so* /mopac/build/
COPY --from=build /mopac/build/mopac /mopac/build/mopac

ENV PATH="$PATH:/mopac/build"

WORKDIR /home
RUN printf '#!/bin/bash\n\
\n\
# Ensure files do not require root \n\
umask 000 \n\
\n\
# Go to mounted volume - /workdir
if [ ! -d /workdir ]; then\n\
  echo "To run mopac, you must mount a volume to /workdir: docker run --rm -v YOURDIR:/workdir ..."\n\
  exit 1\n\
fi\n\
\n\
cd /workdir\n\
"$@"\n'\
>> /home/entry.sh && chmod +x /home/entry.sh

ENTRYPOINT ["/home/entry.sh"]



