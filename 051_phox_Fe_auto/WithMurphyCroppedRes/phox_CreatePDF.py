# Objectives:
# Collect results in tex tables
# Collect figures


# Front matter
##############
import os
from os import fdopen, remove
from tempfile import mkstemp
from shutil import move
import glob
import re
import pandas as pd


# Input information
###################

# PDF file name
pdf_filename = 'Figs_phox_Fe'
file_title = 'phox Results for Fe with auto fitting'
author = 'Rachel Morrison'
composition = 'Fe'


# Functions
###########

# Create header for pdf
def tex_title_info(file_title,author):
	tex_title = '''\\title{'''+file_title+'''}
\\author{'''+author+'''}
\\date	{\\today}'''
	return tex_title

# Create table in LaTeX from .csv file
def latex_table_from_csv(csv_file):
	table_tex = '''\\DTLloadrawdb{db}{'''+csv_file+'''}
\\begin{table}
  \\centering
  \\DTLdisplaydb{db}
\\end{table}
'''
	return table_tex

# Add figures to LaTeX document
def latex_figures(folder,index):
	fig_tex = '''\\begin{figure}[h!]
\\centering
\\subfloat{\\includegraphics[width=8.0cm]{'''+folder+'''/PhoxFigures/NRIXS_'''+index+'''.pdf}} \\quad
\\subfloat{\\includegraphics[width=7.5cm]{'''+folder+'''/PhoxFigures/PeakSub_'''+index+'''.pdf}} \\\\
\\subfloat{\\includegraphics[width=7.5cm]{'''+folder+'''/PhoxFigures/PDOS_'''+index+'''.pdf}} \\quad
\\subfloat{\\includegraphics[width=7.5cm]{'''+folder+'''/PhoxFigures/IntPDOS_'''+index+'''.pdf}} \\\\
\\subfloat{\\includegraphics[width=9.0cm]{'''+folder+'''/PhoxFigures/IntPDOSZoom_'''+index+'''.pdf}}
\\caption{PHOENIX phox data for '''+composition+''' data point '''+index.replace('_','\\_')+'''.}
\\label{fig:cont}
\\end{figure}

'''
	return fig_tex

# Create LaTeX document
def latex_template(title,author,tex_doc):
	template = '''\\documentclass[12pt, oneside]{article}
\\usepackage[top = 1in, bottom = 1in, left = 1in, right = 1in]{geometry}
\\usepackage{titling}
\\usepackage{graphicx}
\\usepackage{subfig}
\\usepackage{datatool}

\\setlength{\\droptitle}{-5em}	% Saves title space, uses package titling

'''+tex_title_info(title,author)+'''

\\begin{document}
\\maketitle
\\tableofcontents

'''+tex_doc+'''

\\end{document}  				
'''
	return template


# Get file structure info
#########################

# Find the filepath of all .res NRIXS files
resfilepath_list = [filepath for filepath in glob.glob('*/*.res')]

folder_list = []
index_dict = dict()
for resfilepath in resfilepath_list:
	folder = re.findall('([A-Za-z0-9/_]+/)[A-Za-z0-9_]+.res',resfilepath)[0]
	index = re.findall('/([A-Za-z0-9_]+).res',resfilepath)[0]
	folder_list.append(folder)
	index_dict[folder] = index


# Create LaTeX Body Text
########################

tex_doc = ''

for folder in folder_list:
	index = index_dict[folder]
	tex_doc = tex_doc+latex_figures(folder,index)


# # Update this part as needed
# folder = '2009Oct_30GPa'
# index = 'Fe_Murphy_P1/'

# tex_doc = latex_figures(folder,index)



# Create pdf
############

build_dir = '.build/'
out_file = build_dir+pdf_filename

tex_code = latex_template(file_title,author,tex_doc)

# Check if build directory exists. Create it if not
if not os.path.exists(build_dir):  # create the build directory if not existing
    os.makedirs(build_dir)

# Create LaTeX document
with open(out_file+'.tex', 'w') as f:  # saves tex_code to output file
    f.write(tex_code)

# Create pdf from LaTeX
os.system('pdflatex -output-directory {} {}'.format(os.path.realpath(build_dir),
	os.path.realpath(out_file)))

move(out_file+'.pdf','Results/'+pdf_filename+'.pdf')