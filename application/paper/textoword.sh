#!/bin/bash

cd statistical_tests_files/figure-latex
for file in *.pdf; do
    mv -- "$file" "${file%%.pdf}"
done
cd ../../

pandoc statistical_tests.tex -o statistical_tests.md
pandoc statistical_tests.md -o statistical_tests.docx
rm statistical_tests.md
