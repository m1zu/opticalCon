TARGET = dummyRootProject
TEMPLATE = app
CONFIG += qt debug warn_on
QMAKE_CXXFLAGS+= -std=c++11

SOURCES += \
  main.cc

CLASSES+= \
  Helpers \
  ProgressBar \
  Functions \
  RootExecFunctionWrapper

for(class, CLASSES) {
  HEADERS += $${class}.hh
  SOURCES += $${class}.cc
}

INCLUDEPATH += \
  $(ROOTSYS)/include

ROOTINCDIR = $$system(root-config --incdir)
!exists ($$ROOTINCDIR/TObject.h) : error("Could NOT find ROOT!")
include ($$ROOTINCDIR/rootcint.pri)
LIBS += \
  $$system(root-config --cflags --libs) \
  -lHistPainter

DESTDIR = builds
OBJECTS_DIR = builds/.tmp
MOC_DIR = builds/.tmp
UI_DIR = builds/.tmp
RCC_DIR = builds/.tmp
