

document :
	Rscript -e "devtools::document()"

test :
	Rscript -e "devtools::test()"

vignette :
	@echo file:///`pwd`/doc
	Rscript -e "devtools::build_vignettes()"

site : document
	echo "pkgdown::build_site()" |R --vanilla

publish : 
	scp -r docs/* logarithmic.net:www/weitrix/
