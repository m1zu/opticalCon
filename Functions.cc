#include "Functions.hh"

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
