
document :
	Rscript -e "devtools::document()"

data :
	Rscript -e 'devtools::load_all("."); source("data-raw/simwei.R")'

check.Renviron :
	curl https://raw.githubusercontent.com/Bioconductor/packagebuilder/master/check.Renviron >check.Renviron

check : document check.Renviron
	rm weitrix_*.tar.gz || echo No old tarball
	time R CMD build .
	R_CHECK_ENVIRON=check.Renviron time R CMD check weitrix_*.tar.gz
	ls -lh weitrix_*.tar.gz

bioccheck : document
	R CMD BiocCheck .

test :
	Rscript -e "devtools::test()"

install :
	R CMD INSTALL .

vignette :
	@echo file:///`pwd`/doc
	Rscript -e "devtools::build_vignettes(quiet=FALSE,clean=FALSE)"

site : document
	Rscript -e "pkgdown::build_site(new_process=FALSE)"

# Note: must install package for system.file to work
site-devel : document
	Rscript -e "devtools::load_all('.',export_all=F);pkgdown::build_site(new_process=FALSE,devel=TRUE,lazy=TRUE)"

# Note: must install package for system.file to work
articles-devel : 
	Rscript -e "devtools::load_all('.',export_all=F);pkgdown::build_articles(lazy=TRUE,quiet=FALSE)"

publish : 
	rsync -rv docs/* doc/* logarithmic.net:www/weitrix/

clean :
	rm -rf weitrix.Rcheck vignettes/*.html vignettes/*_cache vignettes/*_files docs weitrix_*.tar.gz


.PHONY : document data check bioccheck test install vignette site site-devel articles-devel publish clean
