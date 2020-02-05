
document :
	Rscript -e "devtools::document()"

data :
	Rscript -e 'devtools::load_all("."); source("data-raw/simwei.R")'

check : document
	R CMD build .
	R CMD check weitrix_*.tar.gz
	rm weitrix_*.tar.gz

bioccheck : document
	R CMD BiocCheck .

test :
	Rscript -e "devtools::test()"

install :
	R CMD INSTALL .

vignette :
	@echo file:///`pwd`/doc
	Rscript -e "devtools::build_vignettes()"

site : document
	Rscript -e "pkgdown::build_site(new_process=FALSE)"

# Note: must install package for system.file to work
site-devel : document
	Rscript -e "devtools::load_all('.',export_all=F);pkgdown::build_site(new_process=FALSE,devel=TRUE,lazy=TRUE)"

# Note: must install package for system.file to work
articles-devel : 
	Rscript -e "devtools::load_all('.',export_all=F);pkgdown::build_articles(lazy=TRUE,quiet=FALSE)"

publish : 
	rsync -rv docs/* logarithmic.net:www/weitrix/

clean :
	rm -rf weitrix.Rcheck vignettes/*.html vignettes/*_cache vignettes/*_files


.PHONY : document data check bioccheck test install vignette site site-devel articles-devel publish clean
