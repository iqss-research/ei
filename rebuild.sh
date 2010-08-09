#!/usr/bin/sh


if [ -x /usr/bin/R ]
then
    R="/usr/bin/R"
fi


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
echo "* copying over files"
cp -Rf ei cran/


# useless file removal

if [ -d ei/inst/stata ]
then
    echo " * removing /inst/stata"
    rm -rf ei/inst/stata
fi

if [ -d ei/inst/stext ]
then
    echo " * removing ei/inst/text"
    rm -rf ei/inst/text
fi

if [ -d ei/inst/CVS ]
then
    echo " * removing ei/inst/CVS"
    rm -rf ei/inst/CVS
fi

if [ -d ei/inst/doc/CVS ]
then
    echo " * removing ei/inst/doc/CVS"
    rm -rf ei/inst/doc/CVS
fi


echo "* building the package"
$R CMD build ei

echo "* checking the package"
$R CMD check ei*.tar.gz


if [ -x ei.tar.gz ]
then
    echo "* removing ei.tar.gz"
    rm -f ei.tar.gz
fi


#
echo "* compressing package"
tar -czf ei.tar.gz ei/

#
echo "* installing package"
$R CMD INSTALL ei.tar.gz
