#include "Helpers.hh"

#include <TApplication.h>
#include <TVirtualFitter.h>

double medianLightYield(TH1D* h)
{
  return TMath::Median(h->GetNbinsX(), &(h->GetArray()[1])); //full range
  //return TMath::Median(56, &(h->GetArray()[261])); //shortRange
}

double meanLightYield(TH1D* h)
{
  double mean=0;
  int n = h->GetNbinsX();
  double* channelLightYield = &(h->GetArray()[1]);

  // full range
  for (int i=0; i<n; ++i)
    if ((((i+1)%64)==0 || (i%64)==0))
      ; else mean+=channelLightYield[i];
  mean /= double(n-16);

  /*  // short range
  int count = 0;
  for (int i=0; i<n; ++i)
    if ((i>260) && (i<315)) {
      mean+=channelLightYield[i];
      count++;
    }
  mean /= double(count-1); */

  return mean;
}

double meanLightYieldError(TH1D* h)
{
  double mean = meanLightYield(h);
  int n = h->GetNbinsX();
  double* channelLightYield = &(h->GetArray()[1]);
  double variance = 0;

  // full range
  for (int i=0; i<n; ++i)
    if (((i+1)%64)==0 || (i%64)==0)
      ; else variance+=TMath::Power(mean-channelLightYield[i], 2);
  variance /= double(n-17);
  double rms= TMath::Sqrt(variance)/TMath::Sqrt(n-16);

  /* // short range
  int count = 0;
  for (int i=0; i<n; ++i)
    if ((i>260) && (i<315)) {
      variance+=TMath::Power(mean-channelLightYield[i], 2);
      count++;
    }
  variance /= double(count-1);
  double rms= TMath::Sqrt(variance)/TMath::Sqrt(double(count-1)); */

  return rms;
}

double medianLightYieldError(TH1D* h)
{
  return meanLightYieldError(h)*1.2533;
}

double MAD(TH1D* h)
{
  using namespace TMath;
  Double_t median = medianLightYield(h);
  Double_t* tempArray = new Double_t[h->GetNbinsX()];
  for (int i=0; i<h->GetNbinsX(); ++i) {
    tempArray[i] = Abs( h->GetArray()[i+1] -median);
  }
  return Median(h->GetNbinsX(), tempArray);
}

std::vector<TF1*> doubleExpFit_Fresnel(TGraphErrors* g1, TGraphErrors* g2)
{
  std::vector<TF1*> func;
  TF1* f1 = new TF1("f1", "[0]*((1-[2])*exp(-x/[1])+[2]*exp(-x/[3]))+ [0]*[4]*((1-[2])*exp(-(242.35*2-x)/[1])+[2]*exp(-(242.35*2-x)/[3]))", 0., 200.);
  TF1* f2 = new TF1("f2", "[0]*((1-[2])*exp(-x/[1])+[2]*exp(-x/[3]))+ [0]*[4]*((1-[2])*exp(-(242.35*2-x)/[1])+[2]*exp(-(242.35*2-x)/[3]))", 0., 200.);
  f1->SetParameters( 20., 50., 0.7, 400., 0.089);
  f1->SetParNames( "I_{0}", "#Lambda_{s}", "f_{l}", "#Lambda_{l}", "r");
  f1->SetParLimits(1, 1., 800.);
  f1->SetParLimits(3, 1., 800.);
  f1->SetParLimits(2, .6, .9);
  f1->FixParameter(4, 0.089);
  f1->SetLineWidth(0.2);
  f1->SetLineColor(kRed);
  g1->Fit("f1", "MEQ N");

  TGraphErrors* grint1 = new TGraphErrors(250.);
  for (int i=0; i<250.; ++i)
    grint1->SetPoint(i,i,0);
  TVirtualFitter::GetFitter()->GetConfidenceIntervals(grint1,0.68);
  grint1->SetLineColor(kRed);
  grint1->Draw("e5  SAME");


  f2->SetParameters( 20., 50., 0.7, 400., 0.089);
  f2->SetParNames( "I_{0}", "#Lambda_{s}", "f_{l}", "#Lambda_{l}", "r");
  f2->SetParLimits(1, 1., 800.);
  f2->SetParLimits(3, 1., 800.);
  f2->SetParLimits(2, .6, .9);
  f2->FixParameter(4, 0.089);
  f2->SetLineWidth(0.2);
  f2->SetLineColor(kBlue);
  g2->Fit("f2", "MEQ N");

  TGraphErrors* grint2 = new TGraphErrors(250.);
  for (int i=0; i<250.; ++i)
    grint2->SetPoint(i,i,0);
  TVirtualFitter::GetFitter()->GetConfidenceIntervals(grint2,0.68);
  grint2->SetLineColor(kBlue);
  grint2->Draw("e5  SAME");

  {
    double y1= grint1->Eval(30.);
    double y2= grint1->Eval(0.);
    double e1= grint1->GetErrorY(30);
    double e2= grint1->GetErrorY(0);
    //qDebug() << y1 << y2 << e1 << e2;
    qDebug() << " Expected loss on 30cm = " << y1/y2 << "pm  " << TMath::Sqrt((1/y1/y1)*e1*e1 + (y1/y2/y2)*(y1/y2/y2)*e2*e2);
  }

  gStyle->SetOptFit(0);
  /*
  TLegend* l = new TLegend(.55,.38,.89,.89);
  l->AddEntry(expFit, TString::Format("#chi^{2}/ndf  (%.2f/%i)", expFit->GetChisquare(), expFit->GetNDF()), "l");
  l->AddEntry((TObject*)0, TString::Format("I_{0}  (%.2f #pm %.2f)pixel", expFit->GetParameter(0), expFit->GetParError(0)), "");
  l->AddEntry((TObject*)0, TString::Format("f_{long}  (%.2f #pm %.2f)", expFit->GetParameter(2), expFit->GetParError(2)), "");
  l->AddEntry((TObject*)0, TString::Format("#lambda_{long}  (%.1f #pm %.1f)cm", expFit->GetParameter(3), expFit->GetParError(3)), "");
  l->AddEntry((TObject*)0, TString::Format("#lambda_{short}  (%.1f #pm %.1f)cm", expFit->GetParameter(1), expFit->GetParError(1)), "");
  l->AddEntry((TObject*)0, TString::Format("r_{Fresnel} [fix]  %.3f", expFit->GetParameter(4)), "");
  l->Draw("SAME");
  */
  func.push_back(f1);
  func.push_back(f2);
  return func;
}

int main(int argc, char** argv) {
  TApplication application("analysis", &argc, argv);
  //TH1::AddDirectory(false);
  Helpers::setRootStyle();

  /* extracting median, error on median, position from files */

  QString filepath = ("/home/dami/master/opticalCon/4TSAACFIM00127/noMirror/");
  QStringList files = QDir(filepath).entryList(QStringList() << "*AZ*.root");
  const int nFiles = files.count();
  qDebug() << endl << "using" << nFiles << " files before cut.. ";

  Double_t* mean = new Double_t[nFiles];
  Double_t* median = new Double_t[nFiles];
  Double_t* meanError = new Double_t[nFiles];
  Double_t* medianError = new Double_t[nFiles];
  Double_t* position = new Double_t[nFiles];

  int count = 0;
  foreach (QString file, files)
  {
    qDebug() << "\n  ..processing "<< file;

    QString fullPath = filepath + file;
    auto hist = Helpers::extractFromFile<TH1D>(fullPath, "cluster signal value 0 HistogramDisplay").front();
    mean[count] = meanLightYield(hist);
    median[count] = medianLightYield(hist);
    position[count] = file.left(3).toDouble();
    meanError[count] = meanLightYieldError(hist);
    medianError[count] = medianLightYieldError(hist);

    count ++;
  }

  /* distance to SiPM and error calculation */

  Double_t* positionError = new Double_t[nFiles];

  Double_t posError = TMath::Sqrt((3*0.1*0.1 + 0.01*0.01 + 0.02*0.02)/12.);
  for (int i=0; i<nFiles; ++i) {
    position[i] += (0.3 + 8.7 - 6.5 + 0.4);
    positionError[i] = posError;
  }

  /* reading Files with optical Connection */

  filepath = ("/home/dami/master/opticalCon/4TSAACFIM00127/conNoMirror/");
  files = QDir(filepath).entryList(QStringList() << "*posA*.root");
  const int nCon = files.count();
  qDebug() << endl << "using" << nCon << "files with optical connection.. ";

  Double_t* mean2 = new Double_t[nCon];
  Double_t* median2 = new Double_t[nCon];
  Double_t* meanError2 = new Double_t[nCon];
  Double_t* medianError2 = new Double_t[nCon];
  Double_t* position2 = new Double_t[nCon];

  count = 0;
  foreach (QString file, files)
  {
    qDebug() << "\n  ..processing "<< file;

    QString fullPath = filepath + file;
    auto hist = Helpers::extractFromFile<TH1D>(fullPath, "cluster signal value 0 HistogramDisplay").front();
    mean2[count] = meanLightYield(hist);
    median2[count] = medianLightYield(hist);
    position2[count] = file.left(3).toDouble();
    meanError2[count] = meanLightYieldError(hist);
    medianError2[count] = medianLightYieldError(hist);

    count ++;
  }

  /* optical Connection files: distance to SiPM and error calculation */

  Double_t* positionError2 = new Double_t[nCon];

  for (int i=0; i<nCon; ++i) {
    position2[i] += (0.3 + 8.7 - 6.5 + 0.4);
    positionError2[i] = posError;
  }

  /* draw and fit */

  TCanvas* c1 = new TCanvas("c1", "c1");
  c1->cd();
  TGraphErrors* g1 = new TGraphErrors(nFiles, position, mean, positionError, meanError);
  TGraphErrors* g1_con = new TGraphErrors(nCon, position2, mean2, positionError2, meanError2);
  g1->SetTitle("Mean light yield");
  g1->GetXaxis()->SetTitle("distance to SiPM [cm]");
  g1->GetXaxis()->SetTitleOffset(1.25);
  g1->GetXaxis()->SetTitleSize(20);
  g1->GetXaxis()->SetLabelSize(20);
  g1->GetYaxis()->SetTitleOffset(1.15);
  g1->GetYaxis()->SetTitleSize(20);
  g1->GetYaxis()->SetLabelSize(20);
  g1->GetYaxis()->SetTitle("mean channel signal value [pixel]");
  g1->GetYaxis()->SetRangeUser(8.,26.);
  g1->GetXaxis()->SetLimits(0.,230.);
  g1->Draw("ape");
  g1_con->Draw("pe same");

  TLegend* l1 = new TLegend(.12,.15,.4,.33);
  l1->AddEntry(g1 ,"fiberMat", "pe");
  l1->AddEntry(g1_con ,"fiberMat with opticalCon", "pe");
  l1->Draw("SAME");

  std::vector<TF1*> mean_function_vec = doubleExpFit_Fresnel(g1, g1_con);

  TCanvas* c2 = new TCanvas("c2", "c2");
  c2->cd();
  TGraphErrors* g2 = new TGraphErrors(nFiles, position, median, positionError, medianError);
  TGraphErrors* g2_con = new TGraphErrors(nCon, position2, median2, positionError2, medianError2);
  g2->SetTitle("Median light yield");
  g2->GetXaxis()->SetTitle("distance to SiPM [cm]");
  g2->GetXaxis()->SetTitleOffset(1.15);
  g2->GetXaxis()->SetTitleSize(20);
  g2->GetXaxis()->SetLabelSize(20);
  g2->GetYaxis()->SetTitleOffset(1.25);
  g2->GetYaxis()->SetTitleSize(20);
  g2->GetYaxis()->SetLabelSize(20);
  g2->GetYaxis()->SetTitle("median channel signal value [pixel]");
  g2->GetYaxis()->SetRangeUser(8., 26.);
  g2->GetXaxis()->SetLimits(0., 230.);
  g2->Draw("ape");
  g2_con->Draw("pe same");

  TLegend* l2 = new TLegend(.12,.15,.4,.33);
  l2->AddEntry(g2 ,"fiberMat", "pe");
  l2->AddEntry(g2_con ,"fiberMat with opticalCon", "pe");
  l2->Draw("SAME");

  std::vector<TF1*> median_function_vec = doubleExpFit_Fresnel(g2, g2_con);

  /* evaluation optical connection */
  const int nSamples = 180;
  Double_t* ratioMean = new Double_t[nSamples];
  Double_t* posRatio = new Double_t[nSamples];

  for (int i=0; i<nSamples; ++i){
      posRatio[i] = double(i+20);
      ratioMean[i] = mean_function_vec[1]->Eval(double(i+20)) /  mean_function_vec[0]->Eval(double(i+20));
  }

  TGraph* gRatioMean = new TGraph(nSamples, posRatio, ratioMean);
  TCanvas* c3 = new TCanvas("c3", "c3");
  c3->cd();
  gRatioMean->Draw("apl");

  application.Run();
  return 0;
}
