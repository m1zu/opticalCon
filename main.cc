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

  // full range (except border channels)
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

double meanClusterWidth(TH1D* h)
{
    double mean=0;
    int n = h->GetNbinsX();
    double* channelClusterWidth = &(h->GetArray()[1]);

    for (int i=0; i<n; ++i)
      if ((((i+1)%64)==0 || (i%64)==0))
        ; else mean+=channelClusterWidth[i];
    mean /= double(n-16);

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

  QString filepath = ("/home/iwanicki/sw/opticalCon/data/4TSAACFIM00387/noMirror/");
  QStringList files = QDir(filepath).entryList(QStringList() << "*.root");
  const int nFiles = files.count();
  qDebug() << endl << "using" << nFiles << " files before cut.. ";

  Double_t* mean = new Double_t[nFiles];
  Double_t* median = new Double_t[nFiles];
  Double_t* meanError = new Double_t[nFiles];
  Double_t* medianError = new Double_t[nFiles];
  Double_t* position = new Double_t[nFiles];
  Double_t* clusterWidth = new Double_t[nFiles];

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

    auto histw = Helpers::extractFromFile<TH1D>(fullPath, "cluster width value 0 HistogramDisplay").front();
    clusterWidth[count] = meanClusterWidth(histw);

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

  filepath = ("/home/iwanicki/sw/opticalCon/data/4TSAACFIM00387/withConnectionAndFluid/");
  files = QDir(filepath).entryList(QStringList() << "*posA*.root");
  const int nCon = files.count();
  qDebug() << endl << "using" << nCon << "files with optical connection.. ";

  Double_t* mean2 = new Double_t[nCon];
  Double_t* median2 = new Double_t[nCon];
  Double_t* meanError2 = new Double_t[nCon];
  Double_t* medianError2 = new Double_t[nCon];
  Double_t* position2 = new Double_t[nCon];
  Double_t* clusterWidth2 = new Double_t[nCon];

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

    auto histw = Helpers::extractFromFile<TH1D>(fullPath, "cluster width value 0 HistogramDisplay").front();
    clusterWidth2[count] = meanClusterWidth(histw);

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
  g1->GetXaxis()->SetTitleOffset(1.00);
  g1->GetXaxis()->SetTitleSize(30);
  g1->GetXaxis()->SetLabelSize(30);
  g1->GetYaxis()->SetTitleOffset(.8);
  g1->GetYaxis()->SetTitleSize(30);
  g1->GetYaxis()->SetLabelSize(30);
  g1->GetYaxis()->SetTitle("mean channel signal value [pixel]");
  g1->GetYaxis()->SetRangeUser(8.,26.);
  g1->GetXaxis()->SetLimits(0.,230.);
  g1->Draw("ape");
  g1_con->Draw("pe same");

  std::vector<TF1*> mean_function_vec = doubleExpFit_Fresnel(g1, g1_con);

  TLegend* l1 = new TLegend(.55,.65,.9,.9);
  l1->AddEntry(mean_function_vec[0] ,"fiberMat before cut", "l");
  l1->AddEntry(mean_function_vec[1] ,"fiberMat with opticalCon", "l");
  l1->Draw("SAME");

  TCanvas* c2 = new TCanvas("c2", "c2");
  c2->cd();
  TGraphErrors* g2 = new TGraphErrors(nFiles, position, median, positionError, medianError);
  TGraphErrors* g2_con = new TGraphErrors(nCon, position2, median2, positionError2, medianError2);
  g2->SetTitle("Median light yield");
  g2->GetXaxis()->SetTitle("distance to SiPM [cm]");
  g2->GetXaxis()->SetTitleOffset(1.00);
  g2->GetXaxis()->SetTitleSize(30);
  g2->GetXaxis()->SetLabelSize(30);
  g2->GetYaxis()->SetTitleOffset(.8);
  g2->GetYaxis()->SetTitleSize(30);
  g2->GetYaxis()->SetLabelSize(30);
  g2->GetYaxis()->SetTitle("median channel signal value [pixel]");
  g2->GetYaxis()->SetRangeUser(8., 26.);
  g2->GetXaxis()->SetLimits(0., 230.);
  g2->Draw("ape");
  g2_con->Draw("pe same");

  std::vector<TF1*> median_function_vec = doubleExpFit_Fresnel(g2, g2_con);

  TLegend* l2 = new TLegend(.55,.65,.9,.9);
  l2->AddEntry(median_function_vec[0] ,"fiberMat before cut", "l");
  l2->AddEntry(median_function_vec[1] ,"fiberMat with opticalCon", "l");
  l2->Draw("SAME");

  /* evaluation optical connection */
  const int nSamples = 150;
  Double_t* ratioMean = new Double_t[nSamples];
  Double_t* ratioMedian = new Double_t[nSamples];
  Double_t* posRatio = new Double_t[nSamples];

  TF1* shortMat_mean = new TF1("shortMat_mean", "[0]*((1-[2])*exp(-x/[1])+[2]*exp(-x/[3]))+ [0]*[4]*((1-[2])*exp(-(237.35*2-x)/[1])+[2]*exp(-(237.35*2-x)/[3]))", 0., 200.);
  shortMat_mean->SetParameters(mean_function_vec[0]->GetParameter(0),
          mean_function_vec[0]->GetParameter(1),
          mean_function_vec[0]->GetParameter(2),
          mean_function_vec[0]->GetParameter(3),
          mean_function_vec[0]->GetParameter(4));

  TF1* shortMat_median = new TF1("shortMat_median", "[0]*((1-[2])*exp(-x/[1])+[2]*exp(-x/[3]))+ [0]*[4]*((1-[2])*exp(-(237.35*2-x)/[1])+[2]*exp(-(237.35*2-x)/[3]))", 0., 200.);
  shortMat_median->SetParameters(median_function_vec[0]->GetParameter(0),
          median_function_vec[0]->GetParameter(1),
          median_function_vec[0]->GetParameter(2),
          median_function_vec[0]->GetParameter(3),
          median_function_vec[0]->GetParameter(4));

  for (int i=0; i<nSamples; ++i){
      posRatio[i] = double(i+40);
      //double weight_mean = mean_function_vec[0]->Eval(double(i+40)) / shortMat_mean->Eval(double(i+40));
      //double weight_median = median_function_vec[0]->Eval(double(i+40)) / shortMat_median->Eval(double(i+40));
      ratioMean[i] = mean_function_vec[1]->Eval(double(i+40)) /  mean_function_vec[0]->Eval(double(i+40)); // * weight_mean
      ratioMedian[i] = median_function_vec[1]->Eval(double(i+40)) /  median_function_vec[0]->Eval(double(i+40));  // * weight_median
  }

  TGraph* gRatioMean = new TGraph(nSamples, posRatio, ratioMean);
  TCanvas* c3 = new TCanvas("c3", "c3");
  c3->cd();
  gRatioMean->SetTitle("Ratio of signals with/without c onnection piece");
  gRatioMean->GetXaxis()->SetTitle("distance to SiPM / cm");
  gRatioMean->GetXaxis()->SetTitleOffset(1.0);
  gRatioMean->GetXaxis()->SetTitleSize(30);
  gRatioMean->GetYaxis()->SetTitleOffset( .8);
  gRatioMean->GetYaxis()->SetTitleSize(30);
  gRatioMean->GetYaxis()->SetTitle("signal value ratio (mean values)");
  gRatioMean->Draw("ap");

  TGraph* gRatioMedian = new TGraph(nSamples, posRatio, ratioMedian);
  TCanvas* c4 = new TCanvas("c4", "c4");
  c4->cd();
  gRatioMedian->SetTitle("Ratio of signals with/without connection piece");
  gRatioMedian->GetXaxis()->SetTitle("distance to SiPM / cm");
  gRatioMedian->GetXaxis()->SetTitleOffset(1.0);
  gRatioMedian->GetXaxis()->SetTitleSize(30);
  gRatioMedian->GetYaxis()->SetTitleOffset(.8);
  gRatioMedian->GetYaxis()->SetTitleSize(30);
  gRatioMedian->GetYaxis()->SetTitle("signal value ratio (median values)");
  gRatioMedian->Draw("ap");

  TGraph* gCluster = new TGraph(nFiles, position, clusterWidth);
  TGraph* gClusterCon = new TGraph(nCon, position2, clusterWidth2);
  TCanvas* c5 = new TCanvas("c5", "c5");
  c5->cd();
  gCluster->SetTitle("Mean cluster width");
  gCluster->GetXaxis()->SetTitle("distance to SiPM / cm");
  gCluster->GetXaxis()->SetTitleOffset(1.0);
  gCluster->GetXaxis()->SetTitleSize(30);
  gCluster->GetYaxis()->SetTitleOffset(.8);
  gCluster->GetYaxis()->SetTitleSize(30);
  gCluster->GetYaxis()->SetTitle("mean cluster width (all channels)");
  gCluster->GetYaxis()->SetRangeUser(2.2, 3.0);
  gCluster->SetMarkerColor(kRed);
  gClusterCon->SetMarkerColor(kBlue);
  gCluster->Draw("ap");
  gClusterCon->Draw("same p");

  TLegend* l3 = new TLegend(.60,.65,.9,.9);
  l3->AddEntry(gCluster ,"fiberMat before cut", "p");
  l3->AddEntry(gClusterCon ,"fiberMat with opticalCon", "p");
  l3->Draw("SAME");


  application.Run();
  return 0;
}
