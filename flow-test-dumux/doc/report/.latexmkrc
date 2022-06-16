$ENV{'TZ'}='America/Bogota';

$pdf_mode = 1;   # tex -> pdf
# $pdf_mode = 2; # tex -> ps -> pdf
# $postscript_mode = $dvi_mode = 0;

$pdflatex = 'lualatex %O %S -shell-escape -interaction=batchmode -draftmode';
$pdflatex = 'lualatex %O %S -shell-escape -interaction=batchmode -draftmode';
$pdflatex = 'lualatex %O %S -shell-escape -interaction=batchmode -synctex=1';

$clean_ext = "aux fdb_latexmk fls lof log";
$pdf_previewer = 'start zathura';