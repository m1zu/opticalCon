#ifndef HELPERS_HH
#define HELPERS_HH

#include "RootExecFunctionWrapper.hh"
#include "Functions.hh"

#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TMath.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <THStack.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TList.h>
#include <TKey.h>

#include <cassert>
#include <vector>

#include <QtGlobal>
#include <QString>
#include <QStringList>
#include <QVector>
#include <QMap>
#include <QList>
#include <QFile>
#include <QFileInfo>
#include <QDir>
#include <QTextStream>
#include <QTime>
#include <QDateTime>
#include <QRegExp>
#include <QDebug>

class Helpers
{
public:
  enum ProjectionAxis {ProjectionX, ProjectionY};

  static void setRootStyle();

  static void drawProjection(RootPadEvent, TVirtualPad*, TH2*, ProjectionAxis, int nBins);
  static void drawAndFitProjection(RootPadEvent, TVirtualPad*, TH2*, ProjectionAxis, int nBins, InitializableFunction*);

  static int rootColor(int i);
  static int rootMarkerFull(int i);
  static int rootMarkerOpen(int i);

  template <class T>
  static std::vector<T*> extractFromFile(const QString& fileName, const QString& nameFilter = ".*") {
    return  extractFromFile<T>(fileName, QRegExp(nameFilter));
  }

  template <class T>
  static std::vector<T*> extractFromFile(const QString& fileName, const QRegExp& nameFilter) {
    TVirtualPad* pad = gPad;
    TFile file(qPrintable(fileName));
    gROOT->cd();
    std::vector<T*> objects = extractFromObject<T>(&file, nameFilter);
    std::vector<T*> objectClones(objects.size());
    for (unsigned int i = 0; i < objects.size(); ++i)
      objectClones[i] = static_cast<T*>(objects[i]->Clone());
    file.Close();
    if (pad)
      pad->cd();
    return objectClones;
  }

  template <class T>
  static std::vector<T*> extractFromObject(TObject* object, const QRegExp& nameFilter = QRegExp()) {
    std::vector<T*> objects;

    if (dynamic_cast<T*>(object) && QString(object->GetName()).contains(nameFilter))
      qDebug() << "Found" << object->ClassName() << object->GetName() << object->GetTitle();

    if (dynamic_cast<T*>(object) && QString(object->GetName()).contains(nameFilter))
      objects.push_back(static_cast<T*>(object));

    if (object->InheritsFrom("TPad")) {
      appendToVector<T>(objects, extractFromObject<T>(static_cast<TPad*>(object)->GetListOfPrimitives(), nameFilter));
    } else if (object->InheritsFrom("THStack")) {
      appendToVector<T>(objects, extractFromObject<T>(static_cast<THStack*>(object)->GetHists(), nameFilter));
    } else if (object->InheritsFrom("TMultiGraph")) {
      appendToVector<T>(objects, extractFromObject<T>(static_cast<TMultiGraph*>(object)->GetListOfGraphs(), nameFilter));
    } else if (object->InheritsFrom("TFile")) {
      appendToVector<T>(objects, extractFromObject<T>(static_cast<TFile*>(object)->GetListOfKeys(), nameFilter));
    } else if (object->InheritsFrom("TH1")) {
      appendToVector<T>(objects, extractFromObject<T>(static_cast<TH1*>(object)->GetListOfFunctions(), nameFilter));
    } else if (object->InheritsFrom("TGraph")) {
      appendToVector<T>(objects, extractFromObject<T>(static_cast<TGraph*>(object)->GetListOfFunctions(), nameFilter));
    } else if (object->InheritsFrom("TKey")) {
      appendToVector<T>(objects, extractFromObject<T>(static_cast<TKey*>(object)->ReadObj(), nameFilter));
    } else if (object->InheritsFrom("TList")) {
      for (int i = 0; i < static_cast<TList*>(object)->GetSize(); ++i)
        appendToVector<T>(objects, extractFromObject<T>(static_cast<TList*>(object)->At(i), nameFilter));
    }

    return objects;
  }

  template <class T>
  static void setFontSize(T* h, double size, double titleOffsetX = -1, double titleOffsetY = -1, int font = -1) {
    h->GetXaxis()->SetLabelSize(size);
    h->GetYaxis()->SetLabelSize(size);
    h->GetXaxis()->SetTitleSize(size);
    h->GetYaxis()->SetTitleSize(size);
    h->GetXaxis()->SetTitleOffset(titleOffsetX < 0 ? 30. / size : titleOffsetX);
    h->GetYaxis()->SetTitleOffset(titleOffsetY < 0 ? 30. / size : titleOffsetY);
    h->GetYaxis()->SetTickLength(gStyle->GetTickLength("Y"));
    h->GetXaxis()->SetTickLength(h->GetYaxis()->GetTickLength() * gStyle->GetCanvasDefW() / gStyle->GetCanvasDefH());
    if (font > 0) {
      h->GetYaxis()->SetTitleFont(font);
      h->GetXaxis()->SetTitleFont(font);
      h->GetXaxis()->SetLabelFont(font);
      h->GetYaxis()->SetLabelFont(font);
    }
  }
  static void setFontSize(TH2* h, double size, double titleOffsetX = -1, double titleOffsetY = -1, int font = -1);

  static void setPalettePosition(TH2* h, double x1 = 0.905, double x2 = 0.915);

  static void scaleHistogramSlices(TH2*, ProjectionAxis = ProjectionY);

  static TGraphErrors* meanGraphFromHistogram(const TH2* histogram, int minEntries=100, double xError=0, bool alongX=true);
  static TGraphErrors* correlateGraphs(const TGraphErrors* g1, const TGraphAsymmErrors* g2);
  static void combineGraphs(TGraphErrors* graph1, const TGraphErrors* graph2);

  static TH1D* projection(ProjectionAxis, const TH2*, const char* name = "_px", int firstBin = 0, int lastBin = -1, Option_t* = "");

  static TLatex* drawLatex(const QString& text, double x, double y, double fontSize = 16, Qt::Alignment = Qt::AlignLeft | Qt::AlignBottom);
  static void redrawFrameBox();
  static std::vector<TPad*> dividePad(TVirtualPad*, const std::vector<double>& colPositions, const std::vector<double>& rowPositions = {});

  static void writePdfFile(const std::vector<TPad*>&, const QString& fileName);
  static void writeRootFile(const std::vector<TPad*>&, const QString& fileName);
private:
  template <class T>
  static void appendToVector(std::vector<T*>& objects, const std::vector<T*>& newObjects) {
    objects.insert(objects.end(), newObjects.begin(), newObjects.end());
  }
};

#endif
