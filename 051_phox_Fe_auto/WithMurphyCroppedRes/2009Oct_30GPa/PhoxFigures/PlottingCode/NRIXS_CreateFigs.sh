# Given a file name, this program creates 5 figures using xmgrace:
#	1. A plot of the NRIXS and resolution raw data
#	2. A plot of the peak subtraction
#	3. A plot of the phonon density of states (PDOS)
#	4. A plot of the PDOS integrated over energy
#	5. A zoomed-in plot of the PDOS integral

# INPUT parameter
FILE_NAME="Fe_Murphy_P1"


# PLOT 1: NRIXS
# -------------

INPUTFILE_DAT="../../$FILE_NAME.dat"
INPUTFILE_RES="../../$FILE_NAME.res"

OUTPUTFILE_NRIXS="../NRIXS_$FILE_NAME.ps"
PDFFILE_NRIXS="../NRIXS_$FILE_NAME.pdf"

BATCHFILE_NRIXS="NRIXS.bfile"

# Change NRIXS batchfile
# 	awk (use the awk program)
#	-v VAR1=$INPUTFILE_DAT (allows shell variable VAR1 to be passed to awk)
# 	'/READ XYDY/ {$3 = VAR1} {print $0}' (Find the line containing READ XYDY and replace
#		the third word on the line with VAR1)
#	Then do a funny dance to properly save the new file
awk -v VAR1=\"$INPUTFILE_RES\" '/READ XY/ {$3 = VAR1} {print $0}' $BATCHFILE_NRIXS > input_file.tmp && mv input_file.tmp $BATCHFILE_NRIXS
awk -v VAR1=\"$INPUTFILE_DAT\" '/READ XYDY/ {$3 = VAR1} {print $0}' $BATCHFILE_NRIXS > input_file.tmp && mv input_file.tmp $BATCHFILE_NRIXS
awk -v VAR1=\"$OUTPUTFILE_NRIXS\" '/PRINT TO/ {$3 = VAR1} {print $0}' $BATCHFILE_NRIXS > input_file.tmp && mv input_file.tmp $BATCHFILE_NRIXS

# Create the postscript file
xmgrace -batch NRIXS.bfile -nosafe -hardcopy

# Convert postscript file to pdf
pstopdf $OUTPUTFILE_NRIXS

# Eliminate whitespace around pdf (pdfcrop input.pdf output.pdf, overrides file)
pdfcrop $PDFFILE_NRIXS $PDFFILE_NRIXS

# Delete postscript file to keep clutter down
rm $OUTPUTFILE_NRIXS


# PLOT 2: Resolution Peak Subtraction
# -----------------------------------

INPUTFILE_PSN="../../Output/$FILE_NAME\_psn.dat"

OUTPUTFILE_PSN="../PeakSub_$FILE_NAME.ps"
PDFFILE_PSN="../PeakSub_$FILE_NAME.pdf"

BATCHFILE_PSN="PeakSub.bfile"

# Change PeakSub batchfile
awk -v VAR1=\"$INPUTFILE_PSN\" '/READ XYDY/ {$3 = VAR1} {print $0}' $BATCHFILE_PSN > input_file.tmp && mv input_file.tmp $BATCHFILE_PSN
awk -v VAR1=\"$OUTPUTFILE_PSN\" '/PRINT TO/ {$3 = VAR1} {print $0}' $BATCHFILE_PSN > input_file.tmp && mv input_file.tmp $BATCHFILE_PSN

# Create the postscript file
xmgrace -batch PeakSub.bfile -nosafe -hardcopy

# Convert postscript file to pdf
pstopdf $OUTPUTFILE_PSN

# Eliminate whitespace around pdf (pdfcrop input.pdf output.pdf, overrides file)
pdfcrop $PDFFILE_PSN $PDFFILE_PSN

# Delete postscript file to keep clutter down
rm $OUTPUTFILE_PSN


# PLOT 3: Phonon Density of States
# -----------------------------------

INPUTFILE_PDOS="../../Output/$FILE_NAME\_dos.dat"

OUTPUTFILE_PDOS="../PDOS_$FILE_NAME.ps"
PDFFILE_PDOS="../PDOS_$FILE_NAME.pdf"

BATCHFILE_PDOS="PDOS.bfile"

# Change PeakSub batchfile
awk -v VAR1=\"$INPUTFILE_PDOS\" '/READ XYDY/ {$3 = VAR1} {print $0}' $BATCHFILE_PDOS > input_file.tmp && mv input_file.tmp $BATCHFILE_PDOS
awk -v VAR1=\"$OUTPUTFILE_PDOS\" '/PRINT TO/ {$3 = VAR1} {print $0}' $BATCHFILE_PDOS > input_file.tmp && mv input_file.tmp $BATCHFILE_PDOS

# Create the postscript file
xmgrace -batch PDOS.bfile -nosafe -hardcopy

# Convert postscript file to pdf
pstopdf $OUTPUTFILE_PDOS

# Eliminate whitespace around pdf (pdfcrop input.pdf output.pdf, overrides file)
pdfcrop $PDFFILE_PDOS $PDFFILE_PDOS

# Delete postscript file to keep clutter down
rm $OUTPUTFILE_PDOS


# PLOT 4: Integrated PDOS wrt Energy
# ----------------------------------

INPUTFILE_INTPDOS="../../Output/$FILE_NAME\_dos.dat"

OUTPUTFILE_INTPDOS="../IntPDOS_$FILE_NAME.ps"
PDFFILE_INTPDOS="../IntPDOS_$FILE_NAME.pdf"

BATCHFILE_INTPDOS="IntPDOS.bfile"

# Change PeakSub batchfile
awk -v VAR1=\"$INPUTFILE_INTPDOS\" '/READ XY/ {$3 = VAR1} {print $0}' $BATCHFILE_INTPDOS > input_file.tmp && mv input_file.tmp $BATCHFILE_INTPDOS
awk -v VAR1=\"$OUTPUTFILE_INTPDOS\" '/PRINT TO/ {$3 = VAR1} {print $0}' $BATCHFILE_INTPDOS > input_file.tmp && mv input_file.tmp $BATCHFILE_INTPDOS

# Create the postscript file
xmgrace -batch IntPDOS.bfile -nosafe -hardcopy

# Convert postscript file to pdf
pstopdf $OUTPUTFILE_INTPDOS

# Eliminate whitespace around pdf (pdfcrop input.pdf output.pdf, overrides file)
pdfcrop $PDFFILE_INTPDOS $PDFFILE_INTPDOS

# Delete postscript file to keep clutter down
rm $OUTPUTFILE_INTPDOS


# PLOT 5: Zoomed integrated PDOS wrt Energy
# -----------------------------------------

INPUTFILE_INTPDOSZOOM="../../Output/$FILE_NAME\_dos.dat"

OUTPUTFILE_INTPDOSZOOM="../IntPDOSZoom_$FILE_NAME.ps"
PDFFILE_INTPDOSZOOM="../IntPDOSZoom_$FILE_NAME.pdf"

BATCHFILE_INTPDOSZOOM="IntPDOSZoom.bfile"

# Change PeakSub batchfile
awk -v VAR1=\"$INPUTFILE_INTPDOSZOOM\" '/READ XY/ {$3 = VAR1} {print $0}' $BATCHFILE_INTPDOSZOOM > input_file.tmp && mv input_file.tmp $BATCHFILE_INTPDOSZOOM
awk -v VAR1=\"$OUTPUTFILE_INTPDOSZOOM\" '/PRINT TO/ {$3 = VAR1} {print $0}' $BATCHFILE_INTPDOSZOOM > input_file.tmp && mv input_file.tmp $BATCHFILE_INTPDOSZOOM

# Create the postscript file
xmgrace -batch IntPDOSZoom.bfile -nosafe -hardcopy

# Convert postscript file to pdf
pstopdf $OUTPUTFILE_INTPDOSZOOM

# Eliminate whitespace around pdf (pdfcrop input.pdf output.pdf, overrides file)
pdfcrop $PDFFILE_INTPDOSZOOM $PDFFILE_INTPDOSZOOM

# Delete postscript file to keep clutter down
rm $OUTPUTFILE_INTPDOSZOOM
