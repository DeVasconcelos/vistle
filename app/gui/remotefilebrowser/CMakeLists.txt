set(CMAKE_AUTOMOC ON)

set(SOURCES
    remotefiledialog.cpp
    filedialoghelper.cpp
    qsidebar.cpp
    abstractfilesystemmodel.cpp
    #localfilesystemmodel.cpp
    remotefilesystemmodel.cpp
    abstractfileinfogatherer.cpp
    remotefileiconprovider.cpp
    remotefileinfogatherer.cpp
    vistlefileinfogatherer.cpp)

set(HEADERS
    remotefiledialog.h
    remotefiledialog_p.h
    filedialoghelper.h
    qsidebar_p.h
    remotefscompleter_p.h
    abstractfilesystemmodel.h
    #localfilesystemmodel.h
    remotefilesystemmodel.h
    remotefilesystemmodel_p.h
    abstractfileinfogatherer.h
    remotefileiconprovider.h
    remotefileinfogatherer.h
    vistlefileinfogatherer.h)

set(FORMS remotefiledialog.ui)

add_library(vistle_remotefiledialog STATIC ${SOURCES} ${HEADERS})

if(VISTLE_USE_QT5)
    qt5_wrap_ui(UI_SRCS ${FORMS})
    target_link_libraries(vistle_remotefiledialog PRIVATE Qt5::Widgets)
else()
    qt6_wrap_ui(UI_SRCS ${FORMS})
    target_link_libraries(vistle_remotefiledialog PRIVATE Qt6::Widgets)
    qt_disable_unicode_defines(vistle_remotefiledialog)
endif()
target_sources(vistle_remotefiledialog PRIVATE ${UI_SRCS})

target_compile_definitions(vistle_remotefiledialog PRIVATE QT_DISABLE_DEPRECATED_BEFORE=0x000000)

target_include_directories(
    vistle_remotefiledialog
    PRIVATE ${CMAKE_CURRENT_BINARY_DIR}
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(vistle_remotefiledialog PRIVATE vistle_util vistle_core)
