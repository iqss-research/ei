echo "creating tarball \"ei.tar.gz\""
tar -czf ei.tar.gz ei

echo "installing ei package"
R CMD INSTALL ei.tar.gz
