# syntax=docker/dockerfile:1
FROM antismash/standalone-lite:8.0.4
ENV version=0.15.0

#we override the base image's entrypoint for appropriate interaction with Nextflow
ENTRYPOINT ["/usr/bin/env"]
SHELL ["/bin/bash", "-c"]
CMD /bin/bash
