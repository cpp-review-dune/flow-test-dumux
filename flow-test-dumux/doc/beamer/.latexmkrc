$ENV{'TZ'}='America/Bogota';

$pdf_mode = 1;

$pdflatex = 'lualatex %O %S -shell-escape -interaction=batchmode -draftmode';
$pdflatex = 'lualatex %O %S -shell-escape -interaction=batchmode -draftmode';
$pdflatex = 'lualatex %O %S -shell-escape -interaction=batchmode -synctex=1';

$clean_ext = "aux bbl blg fdb_latexmk fls log";
$pdf_previewer = 'start zathura';