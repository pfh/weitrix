

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

site-devel : 
	Rscript -e "pkgdown::build_site(new_process=FALSE,devel=TRUE,install=TRUE,lazy=TRUE)"

articles-devel : 
	Rscript -e "pkgdown::build_articles(new_process=FALSE,devel=TRUE,lazy=TRUE,quiet=FALSE)"


publish : 
	rsync -rv docs/* logarithmic.net:www/weitrix/
