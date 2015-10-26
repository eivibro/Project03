TEMPLATE = app
CONFIG += console c++14
CONFIG -= app_bundle
#CONFIG -= qt

SOURCES += main.cpp \
    lib.cpp \
    gauss-laguerre.cpp

HEADERS += \
    lib.h

QMAKE_CXXFLAGS += -fopenmp
LIBS += -fopenmp
