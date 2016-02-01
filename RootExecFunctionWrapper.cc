#include "RootExecFunctionWrapper.hh"

#include <TPad.h>
#include <TList.h>

#include <cmath>

RootPadEvent::RootPadEvent(double _x, double _y, int _type)
  : type(_type)
  , x(_x)
  , y(_y)
{
}


RootExecFunctionWrapper::RootExecFunctionWrapper(TPad* pad, const std::function<void(RootPadEvent)>& f)
  : TExec("", "")
  , m_callback(f)
{
  // TPad doesn't have the GetListOfExecs() by default. We need to add a dummy
  // and remove it directly. This will create the TPad::fExecs TList.
  if (!pad->GetListOfExecs()) {
    pad->AddExec("", "");
    pad->GetListOfExecs()->RemoveLast();
  }
  pad->GetListOfExecs()->Add(this);
}

void RootExecFunctionWrapper::Exec(const char*) {
  RootPadEvent event;
  event.type = gPad->GetEvent();
  if (!event.type)
    return;
  event.px = gPad->GetEventX();
  event.py = gPad->GetEventY();
  event.x = gPad->AbsPixeltoX(event.px);
  event.y = gPad->AbsPixeltoY(event.py);
  if (gPad->GetLogx())
    event.x = pow(10, event.x);
  if (gPad->GetLogy())
    event.y = pow(10, event.y);
  event.insideFrame = gPad->GetUxmin() < event.x && event.x < gPad->GetUxmax()
    && gPad->GetUymin() < event.y && event.y < gPad->GetUymax();
  m_callback(event);
}
