
# Run unit tests
test:
	R --silent -e 'devtools::test()'

# Run all examples from the docs
run-examples:
	R --silent -e 'devtools::run_examples()'

# Generate test coverage report
coverage:
	R --silent -e 'devtools::test_coverage()'

# Run additional checks for the package
check-package:
	R --silent -e 'devtools::check()'

# Deploy to CRAN
cran-release:
	R --silent -e 'devtools::release()'

# Install a local development package
install:
	R --silent --slave --no-save --no-restore -e 'devtools::install()'

# Build the package and the manual
build:
	R -e 'devtools::build()'
	R -e 'devtools::build_manual()'

# Setup development environment
dev:
	R --silent --slave --no-save --no-restore -e 'install.packages('devtools')'
	R --silent --slave --no-save --no-restore -e 'devtools::install_deps()'
	R --silent --slave --no-save --no-restore -e 'devtools::install_dev_deps()'
