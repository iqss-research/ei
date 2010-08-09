#!/usr/bin/sh

echo "* checking if cran folder exists"

# check if folder exists
if [ -d cran ]
then
    echo "  * cran folder already exists... deleting"
    rm -fr cran
else
    echo "  * folder doesn't exist... moving on"
fi


# copy over
echo " * copying over files"
cp -Rfv cem cran/


# useless file removal

if [ -d cem/inst/stata ]
then
    echo " * removing cem/inst/stata"
    rm -rf cem/inst/stata
fi

if [ -d cem/inst/stext ]
then
    echo " * removing cem/inst/text"
    rm -rf cem/inst/text
fi

if [ -d cem/inst/CVS ]
then
    echo " * removing cem/inst/CVS"
    rm -rf cem/inst/CVS
fi

if [ -d cem/inst/doc/CVS ]
then
    echo " * removing cem/inst/doc/CVS"
    rm -rf cem/inst/doc/CVS
fi


echo "* building the package"
R CMD build cem


echo "* checking the package"
R CMD CHECK cem*.tar.gz






