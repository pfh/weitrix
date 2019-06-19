

document :
	Rscript -e "devtools::document()"

check : document
	R CMD build .
	R CMD check weitrix_*.tar.gz
	rm weitrix_*.tar.gz

test :
	Rscript -e "devtools::test()"

vignette :
	@echo file:///`pwd`/doc
	Rscript -e "devtools::build_vignettes()"

site : document
	Rscript -e "pkgdown::build_site(new_process=FALSE)"

publish : 
	scp -r docs/* logarithmic.net:www/weitrix/
