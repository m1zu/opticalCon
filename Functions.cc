#include "Functions.h"

#include <TH1.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>

#include <QtGlobal>
#include <QDebug>
#include <QString>

#include <cmath>

InitializableFunction::InitializableFunction()
  : m_tf1()
{
}

InitializableFunction::~InitializableFunction()
{
  for (auto x: m_tf1)
    delete x;
}

void InitializableFunction::clear()
{
  for (auto x: m_latexs)
    delete x;
  m_latexs.clear();
  for (auto x: m_legends)
    delete x;
  m_legends.clear();
  for (auto x: m_lines)
    delete x;
  m_lines.clear();
}

void InitializableFunction::draw(Option_t* opt)
{
  for (auto x: m_latexs)
    x->Draw();
  for (auto x: m_legends)
    x->Draw();
  for (auto x: m_lines)
    x->Draw();
  for (auto x: m_tf1)
    x->Draw(opt);
}

// ----

GaussFunction::GaussFunction(const char* name)
  : InitializableFunction()
{
  m_tf1.push_back(new TF1(name, "gaus(0)", -1, 1));
}

bool GaussFunction::fit(TH1* h, double rangeMin, double rangeMax)
{
  clear();
  m_tf1.at(0)->SetRange(rangeMin, rangeMax);
  if (h->GetEntries() > 20) {
    m_tf1.at(0)->SetParameters(h->GetMaximum(), h->GetMean(), h->GetRMS());
    h->Fit(m_tf1.at(0), "QN0R");
    TLatex* latex = new TLatex(.12, .88, qPrintable(QString("#mu = %1").arg(m_tf1.at(0)->GetParameter(1))));
    latex->SetTextAlign(11);
    latex->SetNDC();
    m_latexs.push_back(latex);
    latex = new TLatex(.12, .88, qPrintable(QString("#sigma = %1").arg(m_tf1.at(0)->GetParameter(2))));
    latex->SetTextAlign(11);
    latex->SetNDC();
    m_latexs.push_back(latex);
    return true;
  }
  m_tf1.at(0)->SetParameters(0, 0, 0);
  return false;
}

// ----

LightYieldFit::LightYieldFit(const char* name)
  : InitializableFunction()
  , m_h(0)
  , m_mpv(-1.)
  , m_mean(-1.)
  , m_median(-1.)
{
  m_tf1.push_back(new TF1(name, "[0] * TMath::GammaDist(x, [1], [2], [3])", 1., 2.));
  m_tf1.at(0)->SetLineColor(kRed);
  m_tf1.at(0)->SetParLimits(0, 0, 1e6);
  m_tf1.at(0)->SetParLimits(1, 1.0, 1e6);
  m_tf1.at(0)->SetParLimits(2, 0, 1e-6);
  m_tf1.at(0)->SetParLimits(3, 1e-6, 1e6);

  //m_tf1.push_back(new TF1(name, "[0] * TMath::LogNormal(x, [1], [2], [3])", 1., 2.));
  //m_tf1.at(0)->SetLineColor(kRed);
  //m_tf1.at(0)->SetParLimits(0, 0, 1e6);
  //m_tf1.at(0)->SetParLimits(1, 1e-6, 1.0);
  //m_tf1.at(0)->SetParLimits(2, 0, 1e-6);
  //m_tf1.at(0)->SetParLimits(3, 1e-6, 1e6);
}

LightYieldFit::~LightYieldFit()
{
}

bool LightYieldFit::fit(TH1* h, double rangeMin, double rangeMax)
{
  clear();
  m_mpv = m_mean = m_median = -1.;
  delete m_h;
  m_h = 0;

  m_tf1.at(0)->SetParameters(0, 2.0, -0.0, 1.0); //gamma
  //m_tf1.at(0)->SetParameters(0, 1.0, -1.0, 1.0); //log normal

  if (h->GetEntries() < 100)
    return false;

  int rebinFactor = 10;
  h->SetLineColor(kGray);
  h->SetLineWidth(2);
  m_h = h->Rebin(rebinFactor, qPrintable(h->GetName() + QString("Rebinned")));
  m_h->Scale(1./rebinFactor);

  int axisFirstBin = h->GetXaxis()->GetFirst();
  int axisLastBin = h->GetXaxis()->GetLast();
  h->GetXaxis()->UnZoom();

  m_mean = h->GetMean();

  const double sigmaPercentage = 0.68268949;
  const double probabilities[] = {.5 - .5 * sigmaPercentage, .5, .5 + .5 * sigmaPercentage};
  double quantiles[3];
  h->GetQuantiles(3, quantiles, probabilities);
  m_median = quantiles[1];
  double tertileL = m_median - quantiles[0];
  double tertileR = quantiles[2] - m_median;

  rangeMin = qMax(m_median - 1.5 * tertileL, qMax(2.5, rangeMin));
  rangeMax = qMin(m_median + 2.0 * tertileR, rangeMax);
  rangeMin = m_h->GetXaxis()->GetBinLowEdge(m_h->GetXaxis()->FindBin(rangeMin));
  rangeMax = m_h->GetXaxis()->GetBinUpEdge(m_h->GetXaxis()->FindBin(rangeMax));

  m_tf1.at(0)->SetRange(rangeMin, rangeMax);

  double maximum = h->GetMaximum();
  m_tf1.at(0)->SetParameters(0, 2.0, 0.0, 1.0); // gamma
  //m_tf1.at(0)->SetParameters(maximum, 0.01, -100.0, 100.0); //log normal
  h->Fit(m_tf1.at(0), "LR QN0");

  h->GetXaxis()->SetRange(axisFirstBin, axisLastBin);

  m_mpv = m_tf1.at(0)->GetMaximumX(rangeMin, rangeMax);

  QString latexTitle = "lognormal parametrization:";
  for (int par = 0; par < m_tf1.at(0)->GetNpar(); ++par)
    latexTitle+= QString(" %1").arg(m_tf1.at(0)->GetParameter(par), 0, 'f', 2);
  TLatex* latex = new TLatex(.88, .88, qPrintable(latexTitle));
  latex->SetTextColor(kRed);
  m_latexs.push_back(latex);

  latex = new TLatex(.88, .84, qPrintable(QString("mpv = %1").arg(m_mpv, 5, 'f', 2)));
  latex->SetTextColor(kBlue + 2);
  m_latexs.push_back(latex);

  latex = new TLatex(.88, .80, qPrintable(QString("mean = %1").arg(m_mean, 5, 'f', 2)));
  latex->SetTextColor(kMagenta + 2);
  m_latexs.push_back(latex);

  latex = new TLatex(.88, .76, qPrintable(QString("median = %1").arg(m_median, 5, 'f', 2)));
  latex->SetTextColor(kCyan + 2);
  m_latexs.push_back(latex);

  for (auto l: m_latexs) {
    l->SetTextAlign(33);
    l->SetTextSize(12);
    l->SetTextFont(83);
    l->SetNDC();
  }

  TLine* line = new TLine(m_mpv, 0, m_mpv, 1.05 * maximum);
  line->SetLineColor(kBlue + 2);
  m_lines.push_back(line);

  line = new TLine(m_mean, 0, m_mean, 1.05 * maximum);
  line->SetLineColor(kMagenta + 2);
  m_lines.push_back(line);

  line = new TLine(m_median, 0, m_median, 1.05 * maximum);
  line->SetLineColor(kCyan + 2);
  m_lines.push_back(line);

  return true;
}

void LightYieldFit::draw(Option_t* option)
{
  InitializableFunction::draw(option);
  if (!m_h)
    return;
  m_h->SetLineColor(kGreen + 2);
  m_h->SetLineWidth(3);
  m_h->Draw("SAME HIST");
}

// ----

ResidualFit::ResidualFit(const char* name)
  : InitializableFunction()
  , m_mean(-1.)
  , m_median(-1.)
  , m_stddev(-1.)
  , m_tertile(-1.)
  , m_sigma(-1.)
  , m_innerSigma(-1.)
  , m_outerSigma(-1.)
{
  m_tf1.push_back(new TF1(qPrintable(QString(name) + "_singleGauss"),
    "gaus(0)", -10., 10.));
  m_tf1.at(0)->SetLineColor(kGray);
  m_tf1.at(0)->SetLineWidth(1);
  m_tf1.push_back(new TF1(qPrintable(QString(name) + "_doubleGauss"), ResidualFit::doubleGauss, -10., 10., 5));
  m_tf1.at(1)->SetParLimits(1, 0.9, 1.0);
  m_tf1.at(1)->SetLineColor(kRed);
}

ResidualFit::~ResidualFit()
{
}

double ResidualFit::doubleGauss(double* x, double* p)
{
  return p[0] * (p[1] * TMath::Gaus(x[0], p[2], p[3]) + (1-p[1]) * TMath::Gaus(x[0], p[2], p[4]));
}

double ResidualFit::effectiveSigma() const
{
  return sqrt(m_innerSigmaPercentage * pow(m_innerSigma, 2) + (1. - m_innerSigmaPercentage) * pow(m_outerSigma, 2));
}

bool ResidualFit::fit(TH1* h, double rangeMin, double rangeMax)
{
  clear();
  m_mean = m_median = -1;
  m_stddev = m_tertile = m_sigma = m_innerSigma = m_outerSigma = -1;
  m_tf1.at(0)->SetParameters(0, 0, 0);
  m_tf1.at(1)->SetParameters(0, 0, 0, 0);

  if (h->GetEntries() < 100)
    return false;

  int axisFirstBin = h->GetXaxis()->GetFirst();
  int axisLastBin = h->GetXaxis()->GetLast();
  h->GetXaxis()->UnZoom();

  double maximum = h->GetMaximum();

  m_mean = h->GetMean();
  m_stddev = h->GetStdDev();

  const double sigmaPercentage = 0.68268949;
  const double probabilities[] = {.5 - .5 * sigmaPercentage, .5, .5 + .5 * sigmaPercentage};
  double quantiles[3];
  h->GetQuantiles(3, quantiles, probabilities);
  m_median = quantiles[1];
  m_tertile = .5 * (quantiles[2] - quantiles[0]);

  m_tf1.at(0)->SetRange(rangeMin, rangeMax);
  m_tf1.at(0)->SetParameters(maximum, m_median, m_tertile);
  h->Fit(m_tf1.at(0), "LR QN0");
  m_sigma = m_tf1.at(0)->GetParameter(2);

  m_tf1.at(1)->SetRange(rangeMin, rangeMax);
  m_tf1.at(1)->SetParameters(maximum, 0.975, m_median, m_tertile, 3 * m_tertile);
  h->Fit(m_tf1.at(1), "LR QN0");
  m_innerSigma = m_tf1.at(1)->GetParameter(3);
  m_outerSigma = m_tf1.at(1)->GetParameter(4);
  m_innerSigmaPercentage = m_tf1.at(1)->GetParameter(1);

  h->GetXaxis()->SetRange(axisFirstBin, axisLastBin);

  m_latexs.push_back(new TLatex(.88, .84, "double Gauss #sigma ="));
  m_latexs.push_back(new TLatex(.88, .81, qPrintable(QString("%1(%2%), %3(%4%)")
    .arg(m_innerSigma, 0, 'f', 3).arg(qRound(100 * m_innerSigmaPercentage))
    .arg(m_outerSigma, 0, 'f', 3).arg(qRound(100 * (1.-m_innerSigmaPercentage))))));
  m_latexs.push_back(new TLatex(.88, .78, qPrintable(QString("--> #sigma_{eff} = %1").arg(effectiveSigma(), 0, 'f', 3))));
  m_latexs.push_back(new TLatex(.88, .74, qPrintable(QString("Gauss #sigma = %1").arg(m_sigma, 0, 'f', 3))));
  m_latexs.push_back(new TLatex(.88, .70, qPrintable(QString("stddev = %1").arg(m_stddev, 0, 'f', 3))));
  m_latexs.push_back(new TLatex(.88, .66, qPrintable(QString("tertile = %1").arg(m_tertile, 0, 'f', 3))));

  for (auto l: m_latexs) {
    l->SetTextAlign(33);
    l->SetTextSize(12);
    l->SetTextFont(83);
    l->SetNDC();
  }

  TLine* line = new TLine(m_mean, 0, m_mean, 1.05 * maximum);
  line->SetLineColor(kMagenta + 2);
  m_lines.push_back(line);

  line = new TLine(m_median, 0, m_median, 1.05 * maximum);
  line->SetLineColor(kCyan + 2);
  m_lines.push_back(line);

  return true;
}

// ----

bool VerticalBeamProfileAnalysis::fit(TH1* h, double, double) {
  clear();
  TAxis* a = h->GetXaxis();
  const int n = a->GetNbins();
  for (m_lBin = 1; m_lBin <= n && h->GetBinContent(m_lBin) < m_threshold; ++m_lBin);
  for (m_rBin = n; m_rBin >  0 && h->GetBinContent(m_rBin) < m_threshold; --m_rBin);
  ++m_lBin;
  --m_rBin;
  m_lines.push_back(new TLine(a->GetXmin(), m_threshold, a->GetXmax(), m_threshold));
  m_lines.push_back(new TLine(a->GetBinLowEdge(m_lBin), 0, a->GetBinLowEdge(m_lBin), h->GetMaximum()));
  m_lines.push_back(new TLine(a->GetBinUpEdge(m_rBin), 0, a->GetBinUpEdge(m_rBin), h->GetMaximum()));
  return true;
}
