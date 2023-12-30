QT       += core gui webenginewidgets

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++17

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    main.cpp \
    mainwindow.cpp

HEADERS += \
    mainwindow.h

FORMS += \
    mainwindow.ui



#INCLUDEPATH += D:/vcpkg/installed/x64-windows/include
#LIBS += -LD:/vcpkg/installed/x64-windows/lib -lrtklib
INCLUDEPATH += $$PWD/include
LIBS += -L$$PWD/lib -lrtklib
#添加winmm.lib,ws2_32.lib
LIBS += -lwinmm -lws2_32

RC_ICONS = Icon\favicon.ico

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
