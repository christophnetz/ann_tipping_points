TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        Source.cpp \
        individual.cpp \
        rnd.cpp

HEADERS += \
    cxxopts.hpp \
    individual.h \
    rnd.hpp \
    rndutils.hpp

INCLUDEPATH += C:\Users\p288427\MinGW
LIBS += "C:\Users\p288427\MinGW\lib\libboost_program_options.a"
