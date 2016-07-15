#include "Helpers.hh"

#include <TApplication.h>

int main(int argc, char** argv) {
  TApplication application("analysis", &argc, argv);
  //TH1::AddDirectory(false);
  Helpers::setRootStyle();

  TCanvas* canvas = new TCanvas("canvas", "canvas");
  canvas->cd();

  application.Run();
  return 0;
}
