ecbuild_add_executable(TARGET ijedi.x
    SOURCES ijedi.cc
    LIBS ijedi ufo ioda saber)

ecbuild_add_executable(TARGET ijedi_coupled.x
    SOURCES ijedi_coupled.cc
    LIBS ijedi ufo ioda saber)
