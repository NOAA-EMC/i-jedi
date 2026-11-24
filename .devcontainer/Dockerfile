# Use jcsda/docker as base image
FROM noaaepic/ubuntu26.04-gcc14-openmpi5.0.10:v2.2-jedi-bundle

RUN useradd -ms /bin/bash builder

# Set working directory
WORKDIR /home/builder

# Switch to the non-root user
USER builder

# Set default command
CMD ["/bin/bash"]