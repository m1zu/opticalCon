#include "Helpers.hh"

#include <TColor.h>
#include <TPaletteAxis.h>
#include <TFrame.h>

void Helpers::setRootStyle()
{
  int fontStyle = 43;
  int fontSize = 30;

  gStyle->SetCanvasDefH(1100);
  gStyle->SetCanvasDefW(1100 * 16 / 10);
  gStyle->SetStripDecimals(false);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetTitleH(0.06);
  gStyle->SetTitleW(0.50);
  gStyle->SetTitleX(0.50);
  gStyle->SetTitleY(0.98);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetStatColor(0);
  gStyle->SetStatBorderSize(1);
  gStyle->SetLegendBorderSize(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(111);
  gStyle->SetLineScalePS(1.0);
  gStyle->SetHistLineWidth(3);
  gStyle->SetTitleFontSize(.1);
  gStyle->SetLabelOffset(0, "XYZ");
  gStyle->SetLabelFont(fontStyle, "XYZ");
  gStyle->SetLabelSize(fontSize, "XYZ");
  gStyle->SetStatFont(fontStyle);
  gStyle->SetStatFontSize(fontSize);
  gStyle->SetTextFont(fontStyle);
  gStyle->SetTitleFont(fontStyle, "XYZ");
  gStyle->SetTitleSize(fontSize, "XYZ");
  gStyle->SetTitleOffset(2.5, "XYZ");
  gStyle->SetTickLength(0.010, "YZ");
  gStyle->SetTickLength(gStyle->GetTickLength("Y") * gStyle->GetCanvasDefW() / gStyle->GetCanvasDefH(), "X");
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetGridColor(kGray);
  //gStyle->SetLineStyleString(11, "41 20");
  gStyle->SetGridStyle(kDotted);
  gStyle->SetMarkerStyle(kFullCircle);
  gStyle->SetPaintTextFormat(".2f");
  gStyle->SetNumberContours(99);

  TColor::InitializeColors();
  double stops[9] = { 0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000};
  double red[9]   = { 0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764};
  double green[9] = { 0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832};
  double blue[9]  = { 0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539};
  double alpha = 1.0;
  TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
}

int Helpers::rootColor(int i)
{
  if (i == 3)
    return 1;
  if (0 <= i%16 && i%16 <= 7)
    return i+2;
  switch (i%16) {
    case  8: return 28;
    case  9: return 30;
    case 10: return 33;
    case 11: return 38;
  }
  if (12 <= i%16 && i%16 <= 15)
    return i+28;
  return 1;
}

int Helpers::rootMarkerFull(int i)
{
  switch (i % 5) {
    case 0:  return 20;
    case 1:  return 21;
    case 2:  return 29;
    case 3:  return 33;
    case 4:  return 34;
  }
  return 0;
}

int Helpers::rootMarkerOpen(int i)
{
  switch (i % 5) {
    case 0:  return 24;
    case 1:  return 25;
    case 2:  return 30;
    case 3:  return 27;
    case 4:  return 28;
  }
  return 0;
}

void Helpers::scaleHistogramSlices(TH2* histogram, ProjectionAxis axis)
{
  int nBins1 = (axis == ProjectionY) ? histogram->GetXaxis()->GetNbins() : histogram->GetYaxis()->GetNbins();
  int nBins2 = (axis == ProjectionY) ? histogram->GetYaxis()->GetNbins() : histogram->GetXaxis()->GetNbins();
  for (int bin1 = 1; bin1 <= nBins1; ++bin1) {
    TH1D* slice = projection(axis, histogram, (axis == ProjectionY) ? "_py" : "_px", bin1, bin1);
    double maximum = slice->GetMaximum();
    if (!qIsNull(maximum)) {
      for (int bin2 = 1; bin2 <= nBins2 + 1; ++bin2) {
        int binX = (axis == ProjectionY) ? bin1 : bin2;
        int binY = (axis == ProjectionY) ? bin2 : bin1;
        histogram->SetBinContent(binX, binY, histogram->GetBinContent(binX, binY) / maximum);
        histogram->SetBinError(  binX, binY, histogram->GetBinError(  binX, binY) / maximum);
        histogram->SetBinContent(binX, 0, maximum);
      }
    }
  }
}

TGraphErrors* Helpers::meanGraphFromHistogram(const TH2* histogram, int minEntries, double xError, bool alongX)
{
  TGraphErrors* graph = new TGraphErrors();
  const TAxis* axis = alongX ? histogram->GetXaxis() : histogram->GetYaxis();
  for (int bin = 1; bin <= axis->GetNbins(); ++bin) {
    TH1D* slice = projection(alongX ? ProjectionX : ProjectionY, histogram, alongX ? "_px" : "_py", bin, bin);
    if (slice->GetEntries() >= minEntries) {
      double x = axis->GetBinCenter(bin);
      double xErr = xError < 0. ? axis->GetBinWidth(bin) : xError;
      int nPoints = graph->GetN();
      graph->SetPoint(     nPoints, x,    slice->GetMean());
      //graph->SetPoint(     nPoints, x,    slice->GetXaxis()->GetBinCenter(slice->GetMaximumBin()));
      graph->SetPointError(nPoints, xErr, slice->GetMeanError());
      //graph->SetPointError(nPoints, xErr, slice->GetXaxis()->GetBinWidth(slice->GetMaximumBin())); //TODO
    }
  }
  graph->SetTitle(TString::Format("mean_%s;%s;mean %s", histogram->GetTitle(), histogram->GetXaxis()->GetTitle(), histogram->GetYaxis()->GetTitle()));
  return graph;
}

void Helpers::setFontSize(TH2* h, double size, double titleOffsetX, double titleOffsetY, int font)
{
  setFontSize<TH2>(h, size, titleOffsetX, titleOffsetY, font);
  h->GetZaxis()->SetLabelSize(size);
  h->GetZaxis()->SetTitleSize(size);
  h->GetZaxis()->SetTitleOffset(titleOffsetY);
  if (font > 0) {
    h->GetZaxis()->SetTitleFont(font);
    h->GetZaxis()->SetLabelFont(font);
  }
  TPaletteAxis* palette = static_cast<TPaletteAxis*>(h->GetListOfFunctions()->FindObject("palette"));
  if (palette) {
    palette->SetLabelSize(size);
    palette->SetTitleSize(size);
    palette->SetTitleOffset(titleOffsetY);
    if (font > 0) {
      palette->GetAxis()->SetLabelFont(font);
      palette->GetAxis()->SetTitleFont(font);
    }
  }
}

void Helpers::setPalettePosition(TH2* h, double x1, double x2)
{
  gPad->Update();
  TPaletteAxis* palette = static_cast<TPaletteAxis*>(h->GetListOfFunctions()->FindObject("palette"));
  if (palette) {
    palette->GetAxis()->SetTickSize(gStyle->GetTickLength("Y"));
    palette->SetX1NDC(x1);
    palette->SetX2NDC(x2);
  }
}

void Helpers::drawProjection(RootPadEvent event, TVirtualPad* destinationPad, TH2* h, Helpers::ProjectionAxis axis, int nBins)
{
  drawAndFitProjection(event, destinationPad, h, axis, nBins, nullptr);
}

void Helpers::drawAndFitProjection(RootPadEvent event, TVirtualPad* destinationPad, TH2* h, Helpers::ProjectionAxis axis, int nBins, InitializableFunction* f)
{
  if (event.type != kMouseMotion || !event.insideFrame)
    return;

  TH1* projection = 0;
  if (axis == ProjectionX) {
    int bin = h->GetYaxis()->FindBin(event.y);
    projection = Helpers::projection(ProjectionX, h, "_px", bin, bin + nBins - 1);
  }
  if (axis == ProjectionY) {
    int bin = h->GetXaxis()->FindBin(event.x);
    projection = Helpers::projection(ProjectionY, h, "_py", bin, bin + nBins - 1);
  }
  assert(projection);
  if (f)
    f->fit(projection);
  TVirtualPad* previousPad = gPad;
  destinationPad->cd();
  setFontSize(projection, 14);
  projection->Draw("HIST");
  if (f)
    f->draw("SAME");
  gPad->Update();
  if (previousPad)
    previousPad->cd();
}

TGraphErrors* Helpers::correlateGraphs(const TGraphErrors* graph1, const TGraphAsymmErrors* graph2)
{
  TGraphErrors* combinedGraph = new TGraphErrors();
  combinedGraph->SetTitle(TString::Format("%s_vs_%s;%s;%s", graph1->GetTitle(), graph2->GetTitle(), graph1->GetXaxis()->GetTitle(), graph2->GetXaxis()->GetTitle()));
  double* x1 = graph1->GetX();
  double* x2 = graph2->GetX();
  for (int point1 = 0; point1 < graph1->GetN(); ++point1) {
    int correspondingPoint = -1;
    for (int point2 = 0; point2 < graph2->GetN(); ++point2) {
      if (qFuzzyCompare(x1[point1], x2[point2])) {
        correspondingPoint = point2;
        break;
      }
    }
    if (correspondingPoint >= 0) {
      combinedGraph->SetPoint(combinedGraph->GetN(), graph1->GetY()[point1], graph2->GetY()[correspondingPoint]);
    }
  }
  return combinedGraph;
}

void Helpers::combineGraphs(TGraphErrors* graph1, const TGraphErrors* graph2)
{
  for (int i = 0; i < graph2->GetN(); ++i) {
    double x, y;
    graph2->GetPoint(i, x, y);
    graph1->SetPoint(graph1->GetN(), x, y);
//    graph1->SetPointError(graph1->GetN(), graph2->GetErrorX(i), graph2->GetErrorY(i));
  }
}

TH1D* Helpers::projection(ProjectionAxis axis, const TH2* h, const char* name, int firstBin, int lastBin, Option_t* option)
{
  TH1D* projection = 0;
  if (axis == ProjectionX) {
    int f = h->GetXaxis()->GetFirst();
    int l = h->GetXaxis()->GetLast();
    const_cast<TAxis*>(h->GetXaxis())->SetRange();
    projection = h->ProjectionX(name, firstBin, lastBin, option);
    projection->SetTitle(qPrintable(QString("bin %1..%2 #rightarrow %3..%4")
      .arg(firstBin).arg(lastBin)
      .arg(h->GetYaxis()->GetBinLowEdge(firstBin), 0, 'f', 2)
      .arg(h->GetYaxis()->GetBinUpEdge(lastBin), 0, 'f', 2)));
    projection->GetXaxis()->SetRange(f, l);
    const_cast<TAxis*>(h->GetXaxis())->SetRange(f, l);
  } else if (axis == ProjectionY) {
    int f = h->GetYaxis()->GetFirst();
    int l = h->GetYaxis()->GetLast();
    const_cast<TAxis*>(h->GetYaxis())->SetRange();
    projection = h->ProjectionY(name, firstBin, lastBin, option);
    projection->SetTitle(qPrintable(QString("bin %1..%2 #rightarrow %3..%4")
      .arg(firstBin).arg(lastBin)
      .arg(h->GetXaxis()->GetBinLowEdge(firstBin), 0, 'f', 2)
      .arg(h->GetXaxis()->GetBinUpEdge(lastBin), 0, 'f', 2)));
    projection->GetXaxis()->SetRange(f, l);
    const_cast<TAxis*>(h->GetYaxis())->SetRange(f, l);
  } else {
    assert(false);
  }
  return projection;
}

void Helpers::redrawFrameBox()
{
  gPad->Modified();
  gPad->Update();
  TFrame* f = gPad->GetFrame();
  if (!f)
    return;
  TBox* box = new TBox(f->GetX1(), f->GetY1(), f->GetX2(), f->GetY2());
  box->SetFillStyle(0);
  box->Draw();
}

std::vector<TPad*> Helpers::dividePad(TVirtualPad* mother, const std::vector<double>& colPositions, const std::vector<double>& rowPositions)
{
  std::vector<TPad*> newPads;
  if (!rowPositions.size() && !colPositions.size())
    return newPads;
  TVirtualPad* previousPad = gPad;
  mother->Clear();
  mother->cd();
  for (unsigned int row = 1; row < rowPositions.size(); ++row)
    assert(rowPositions[row - 1] < rowPositions[row]);
  for (unsigned int col = 1; col < colPositions.size(); ++col)
    assert(colPositions[col - 1] < colPositions[col]);
  int i = 0;
  for (unsigned int row = 0; row <= rowPositions.size(); ++row) {
    for (unsigned int col = 0; col <= colPositions.size(); ++col) {
      ++i;
      QString name = QString("%1_%2").arg(mother->GetName()).arg(i);
      double x1 = (col == 0) ? 0 : colPositions[col - 1];
      double y1 = (row == rowPositions.size()) ? 0 : 1. - rowPositions[row];
      double x2 = (col == colPositions.size()) ? 1 : colPositions[col];
      double y2 = (row == 0) ? 1 : 1. - rowPositions[row - 1];
      assert(0. <= x1 && x1 <= 1.);
      assert(0. <= x2 && x2 <= 1.);
      assert(0. <= y1 && y1 <= 1.);
      assert(0. <= y2 && y2 <= 1.);
      assert(x1 < x2);
      assert(y1 < y2);
      TPad* pad = new TPad(qPrintable(name), qPrintable(name), x1, y1, x2, y2);
      pad->SetNumber(i);
      pad->Draw();
      newPads.push_back(pad);
    }
  }
  previousPad->cd();
  return newPads;
}

TLatex* Helpers::drawLatex(const QString& text, double x, double y, double fontSize, Qt::Alignment alignment)
{
  TLatex* latex = new TLatex(x, y, qPrintable(text));
  latex->SetNDC();
  latex->SetTextFont(83);
  int h = 0;
  if (alignment.testFlag(Qt::AlignLeft)) h = 1;
  if (alignment.testFlag(Qt::AlignHCenter)) h = 2;
  if (alignment.testFlag(Qt::AlignRight)) h = 3;
  int v = 0;
  if (alignment.testFlag(Qt::AlignTop)) v = 1;
  if (alignment.testFlag(Qt::AlignVCenter)) v = 2;
  if (alignment.testFlag(Qt::AlignBottom)) v = 3;
  latex->SetTextAlign(10 * h + v);
  latex->SetTextSize(fontSize);
  latex->Draw();
  return latex;
}

void Helpers::writePdfFile(const std::vector<TPad*>& canvases, const QString& fileName)
{
  for (unsigned int i = 0; i < canvases.size(); ++i) {
    QString option = fileName;
    if (canvases.size() > 1) {
      if (i == 0) {
        option+= "(";
      } else if (i == canvases.size() - 1) {
        option+= ")";
      }
    }
    canvases[i]->Print(qPrintable(option));
  }
}

void Helpers::writeRootFile(const std::vector<TPad*>& canvases, const QString& fileName)
{
  TFile file(qPrintable(fileName), "RECREATE");
  for (auto canvas: canvases)
    file.WriteTObject(canvas);
  file.Close();
}
