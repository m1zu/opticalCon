#include "Helpers.hh"

#include <TApplication.h>
#include <QDebug>

int main(int argc, char** argv) {
  TApplication application("analysis", &argc, argv);
  Helpers::setRootStyle();

  TCanvas* canvas = new TCanvas("canvas", "canvas");
  canvas->cd();

  application.Run();
  return 0;
}
