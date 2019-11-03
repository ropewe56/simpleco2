#!/bin/bash

pdflatex -output-directory=aux radtrans
biber radtrans
pdflatex -output-directory=aux radtrans

