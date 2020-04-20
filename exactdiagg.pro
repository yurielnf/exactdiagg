TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -larmadillo -lfftw3

SOURCES += main.cpp \
    examples/example_tb.cpp \
    examples/example_hubbard.cpp

HEADERS += \
    exactdiagg/hamsolver.h \
    exactdiagg/qoperator.h \
    exactdiagg/fockbasis.h \
    exactdiagg/symmetrygroup.h \
    exactdiagg/all.h \
    examples/all.h

DISTFILES += \
    examples/plot_hubbard.gnuplot \
    README.md
