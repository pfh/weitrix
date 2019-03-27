

document :
	Rscript -e "devtools::document()"

check : document
	R CMD build .
	R_CHECK_ENVIRON=check.Renviron R CMD check weitrix_*.tar.gz
	rm weitrix_*.tar.gz

test :
	Rscript -e "devtools::test()"

vignette :
	@echo file:///`pwd`/doc
	Rscript -e "devtools::build_vignettes()"

site : document
	echo "pkgdown::build_site()" |R --vanilla

publish : 
	scp -r docs/* logarithmic.net:www/weitrix/
