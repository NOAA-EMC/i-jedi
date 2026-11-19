/*
 * Coupled application: fv3-jedi (atmosphere) + ijedi-mom6 (ocean).
 *
 * This is a TEMPORARY driver used while the ijedi fv3 interface is developed;
 * fv3-jedi will be deprecated here once that interface is ready. It is built
 * only when fv3-jedi is available (BUILD_FV3JEDI=ON in the bundle) and is kept
 * entirely separate from the pure-ijedi mains so that default builds, which do
 * not depend on fv3-jedi, are unaffected.
 */

#include <functional>
#include <map>
#include <memory>
#include <set>
#include <string>

#include "ijedi/Traits.h"

#include "fv3jedi/Utilities/Traits.h"

#include "oops/coupled/GetValuesCoupled.h"
#include "oops/coupled/TraitCoupled.h"

#include "oops/runs/HofX3D.h"
#include "oops/runs/Run.h"
#include "oops/runs/Variational.h"

#include "saber/coupled/instantiateCoupledCovarFactory.h"
#include "saber/oops/instantiateCovarFactory.h"

#include "ufo/instantiateObsFilterFactory.h"
#include "ufo/ObsTraits.h"

// -------------------------------------------------------------------------------------------------

// Atmosphere = real fv3-jedi; ocean = ijedi mom6 (TraitsOcn).
typedef oops::TraitCoupled<fv3jedi::Traits, ijedi::TraitsOcn> CoupledTraits;

// -------------------------------------------------------------------------------------------------

int runApp(int argc, char **argv, const std::string appName)
{
  // Create the Run object
  oops::Run run(argc, argv);

  // Instantiate factories.
  ufo::instantiateObsFilterFactory();
  // The coupled background error covariance is block-per-component. This registers
  // the coupled covariance models ("SABER coupled", "Coupled Block Diagonal") as
  // well as the per-component SABER factories they delegate to.
  saber::instantiateCoupledCovarFactory<fv3jedi::Traits, ijedi::TraitsOcn>();

  // Map from app names to factory lambdas
  std::map<std::string, std::function<std::unique_ptr<oops::Application>()>> apps;

  apps["hofx3d"] = []()
  {
    return std::make_unique<oops::HofX3D<CoupledTraits, ufo::ObsTraits>>();
  };

  apps["var"] = []()
  {
    return std::make_unique<oops::Variational<CoupledTraits, ufo::ObsTraits>>();
  };

  // Create and run the requested application
  auto it = apps.find(appName);
  return run.execute(*(it->second()));
}

// -------------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
  // Check that the number of arguments is correct
  ASSERT_MSG(argc >= 2, "Usage: " + std::string(argv[0]) + " <app> <options>");

  // Get the application to be run
  std::string appName = argv[1];
  for (char &c : appName)
  {
    c = std::tolower(c);
  }

  // Check that the application is recognized
  const std::set<std::string> validApps = {
      "hofx3d",
      "var",
  };
  ASSERT_MSG(validApps.find(appName) != validApps.end(), "Application not recognized: " + appName);

  // Remove program from argc and argv
  argv[1] = argv[0];  // Move executable name to second position
  argv += 1;          // Move pointer up one
  argc -= 1;          // Remove 1 from count

  // Call application specific main functions
  return runApp(argc, argv, appName);
}

// -------------------------------------------------------------------------------------------------
