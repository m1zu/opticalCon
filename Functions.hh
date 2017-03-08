#ifndef InitializableFunction_HH
#define InitializableFunction_HH

#include <TF1.h>

#include <climits>
#include <cfloat>
#include <vector>

class TH1;
class TLatex;
class TLegend;
class TLine;

class InitializableFunction
{
public:
  InitializableFunction();
  virtual ~InitializableFunction();
  virtual bool fit(TH1*, double rangeMin = -DBL_MAX, double rangeMax = DBL_MAX) = 0;
  virtual void draw(Option_t*);
  int numberOfTf1() const {return m_tf1.size();}
  TF1* tf1(int i = 0) const {return m_tf1.at(i);}
protected:
  void clear();
  std::vector<TF1*> m_tf1;
  std::vector<TLatex*> m_latexs;
  std::vector<TLegend*> m_legends;
  std::vector<TLine*> m_lines;
};

class GaussFunction : public InitializableFunction
{
public:
  GaussFunction(const char* name);
  virtual bool fit(TH1*, double rangeMin = -DBL_MAX, double rangeMax = DBL_MAX) override;
};

#endif
