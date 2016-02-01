#ifndef ROOTEXECFUNCTIONWRAPPER_HH
#define ROOTEXECFUNCTIONWRAPPER_HH

#include <TExec.h>
#include <Buttons.h>

#include <functional>

class TPad;

struct RootPadEvent
{
  RootPadEvent(double x = 1, double y = 1, int type = kMouseMotion);
  int type = kMouseMotion;
  int px = 0, py = 0;
  double x = 1, y = 1;
  bool insideFrame = true;
};

class RootExecFunctionWrapper : public TExec
{
public:
  template<class FunctionType, typename... FunctionArgumentTypes>
  static RootExecFunctionWrapper* Create(TPad* pad, FunctionType&& function, FunctionArgumentTypes&&... args) {
    return new RootExecFunctionWrapper(pad, std::bind(function, std::placeholders::_1, args...));
  }
  void trigger(const RootPadEvent event) {m_callback(event);}
private:
  RootExecFunctionWrapper(TPad* pad, const std::function<void(RootPadEvent)>& f);
  virtual void Exec(const char*);
  std::function<void(const RootPadEvent)> m_callback;
};

#endif
