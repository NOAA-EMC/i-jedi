ecbuild_add_executable(TARGET ijedi.x
    SOURCES ijedi.cc
    LIBS ijedi ufo ioda saber)

ecbuild_add_executable(TARGET ijedi_coupled.x
    SOURCES ijedi_coupled.cc
    LIBS ijedi ufo ioda saber)

# Temporary fv3-jedi (atmosphere) + ijedi-mom6 (ocean) coupled driver.
# Built only when fv3-jedi is available (bundle option BUILD_FV3JEDI=ON), so
# default builds without fv3-jedi are unaffected.
if(fv3jedi_FOUND)
  ecbuild_add_executable(TARGET ijedi_fv3_mom6_coupled.x
      SOURCES ijedi_fv3_mom6_coupled.cc
      LIBS ijedi fv3jedi ufo ioda saber)
endif()
