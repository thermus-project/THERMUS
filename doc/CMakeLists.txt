# Cmake for latex documentation
INCLUDE(${CMAKE_SOURCE_DIR}/cmake/modules/UseLATEX.cmake)

ADD_LATEX_DOCUMENT("thermus_reference_v4.tex" 
    BIBFILES "thermus_biblio_v4.bib" 
    IMAGE_DIRS figures 
    DEPENDS elsarticle.cls FORCE_PDF)

#ADD_LATEX_DOCUMENT("testdoc.tex")
