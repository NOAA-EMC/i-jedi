#include <functional>
#include <map>

#include "ijedi/Traits.h"

#include "oops/runs/Run.h"
#include "oops/runs/HofX3D.h"
#include "oops/runs/Variational.h"
#include "oops/runs/ConvertState.h"

#include "saber/oops/instantiateCovarFactory.h"
#include "saber/oops/ErrorCovarianceToolbox.h"

#include "ufo/instantiateObsFilterFactory.h"
#include "ufo/ObsTraits.h"

// -------------------------------------------------------------------------------------------------

int runApp(int argc, char **argv, const std::string appName)
{
  // Create the Run object
  oops::Run run(argc, argv);

  // Test application pointer
  std::unique_ptr<oops::Application> app;

  // Intantiate factories
  ufo::instantiateObsFilterFactory();
  saber::instantiateCovarFactory<ijedi::Traits>();

  // Define a map from app names to lambda functions that create unique_ptr to Applications
  std::map<std::string, std::function<std::unique_ptr<oops::Application>()>> apps;

  apps["hofx3d"] = []() {
    return std::make_unique<oops::HofX3D<ijedi::Traits, ufo::ObsTraits>>();
  };
  apps["var"] = []() {
    return std::make_unique<oops::Variational<ijedi::Traits, ufo::ObsTraits>>();
  };
  apps["convertstate"] = []() {
    return std::make_unique<oops::ConvertState<ijedi::Traits>>();
  };
  apps["errortoolbox"] = []() {
    return std::make_unique<saber::ErrorCovarianceToolbox<ijedi::Traits>>();
  };

  // Create application object and point to it
  auto it = apps.find(appName);

  // Run the application
  return run.execute(*(it->second()));
}

// -------------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
  // Check that the number of arguments is correct
  // ----------------------------------------------
  ASSERT_MSG(argc >= 2, "Usage: " + std::string(argv[0]) + " <app> <options>");

  // Get the application to be run
  std::string appName = argv[1];
  for (char &c : appName)
  {
    c = std::tolower(c);
  }

  // Check that the application is recognized
  // ----------------------------------------
  const std::set<std::string> validApps = {
      "hofx3d", "var", "errortoolbox", "convertstate"
  };
  ASSERT_MSG(validApps.find(appName) != validApps.end(), "Application not recognized: " + appName);

  // Remove program from argc and argv
  // ---------------------------------
  argv[1] = argv[0];  // Move executable name to second position
  argv += 1;          // Move pointer up one
  argc -= 1;          // Remove 1 from count

  // Call application specific main functions
  // ----------------------------------------
  return runApp(argc, argv, appName);
}

// -------------------------------------------------------------------------------------------------
