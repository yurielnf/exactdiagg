TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -larmadillo -lfftw3

SOURCES += main.cpp \
    examples/example_hernan/example_hernan.cpp \
    examples/example_hernan/parameters.cpp \
    examples/example_tb.cpp \
    examples/example_hubbard.cpp \
    examples/example_hubbard2.cpp

HEADERS += \
    exactdiagg/hamsolver.h \
    exactdiagg/qoperator.h \
    exactdiagg/fockbasis.h \
    exactdiagg/symmetrygroup.h \
    exactdiagg/all.h \
    examples/all.h \
    examples/example_hernan/parameters.h

DISTFILES += \
    examples/plot_hubbard.gnuplot \
    README.md
